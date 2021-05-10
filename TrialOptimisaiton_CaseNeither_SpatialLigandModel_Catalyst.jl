using Random, Distributions
using CSV, DataFrames, DelimitedFiles
using Catalyst
using Profile

function wrappeddist(X1,X2,N)

  if abs(X1-X2) >= ceil(Int,N/2)
      Dist = N - abs(X1-X2)
  else
      Dist = abs(X1-X2)
  end

  Dist 
end

function calculate_propensity(State,SpeciesMap,Reactions,ParameterMap,Pos,ParticleType,ReactionType,Nhbd,DoiRateGridCentre,DoiRateMatrix,Doicutoffd,GaussRateGridCentre,GaussRateMatrix,Gausscutoffd)


  Propensities = zeros(Float64,cumsum_RN[end])
  # add the two diffusive reactions
  Propensities[1] = 4*Dreceptor/h^2 * sum(State[ReceptorSpecies])
  Propensities[2] = 4*Dligand/h^2 * sum(State[LigandSpecies])

  for (rx,reac) in enumerate(Reactions)
    if ReactionType[rx+2] == "First order place none" || ReactionType[rx+2] == "First order place one" || ReactionType[rx+2] == "First order place two"
      # we calculate the propensity rates for a first order reactions
      for (tx,t) in enumerate(ParticleType)
        if t==SpeciesMap[reac.substrates[1].op]
          Propensities[cumsum_RN[rx+1]+tx] = ParameterMap[reac.rate.op]
        end
      end
    elseif ReactionType[rx+2] == "Bimolecular Doi" || ReactionType[rx+2] == "Bimolecular Doi Ligand"
      # we calculate the propensity rates for Bimolecular Doi reactions
      # find all particles with type of the first substrate
      T = findall(x->x==SpeciesMap[reac.substrates[1].op], ParticleType)
      for (tx,t) in enumerate(T)
        # find all the particles that are within the Doi reach
        S = findall(x->x==2,Nhbd[t,:])
          for (sx,s) in enumerate(S)
            # check to see if the molecule within reach is of the type of the 2nd substrate
            if ParticleType[s]==SpeciesMap[reac.substrates[2].op] 
              xsep = wrappeddist(Pos[t,1],Pos[s,1],N)
              ysep = wrappeddist(Pos[t,2],Pos[s,2],N)
              Propensities[cumsum_RN[rx+1]+(t-1)*NTotal+s] = ParameterMap[reac.rate.op] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
            end
          end
      end

    elseif ReactionType[rx+2] == "Bimolecular Gauss"
        # we calculate the propensity rates for Bimolecular Doi reactions
      # find all particles with type of the first substrate
      T = findall(x->x==SpeciesMap[reac.substrates[1].op], ParticleType)
      for (tx,t) in enumerate(T)
        # find all the particles that are within the Gassian reach
        S = findall(x->x==1,Nhbd[t,:])
          for (sx,s) in enumerate(S)
            # check to see if the molecule within reach is of the type of the 2nd substrate
            if ParticleType[s]==SpeciesMap[reac.substrates[2].op] 
              xsep = wrappeddist(Pos[t,1],Pos[s,1],N)
              ysep = wrappeddist(Pos[t,2],Pos[s,2],N)
              Propensities[cumsum_RN[rx+1]+(t-1)*NTotal+s] = ParameterMap[reac.rate.op] * GaussRateMatrix[GaussRateGridCentre+xsep,GaussRateGridCentre+ysep]
            end
          end
      end
    end
  end
  return sum(Propensities),Propensities
end

function sample_reaction(a0,Propensities)

  # calculate the cumulative sum as you go
  reacted=0
  cumrate=0
  RX=1
  r = rand()
  while reacted==0
      cumrate += Propensities[RX]
      if a0*r < cumrate
          reacted=1
      else
      RX+=1
      end
  end

  return RX

end

function findRec_or_Lig(ParticleType,ReceptorSpecies)
  Receptors = []
  for (rx,r) in enumerate(ReceptorSpecies)
    if r==1
      X = findall(x->x==rx,ParticleType)
      if length(X)>0
        Receptors = vcat(Receptors,X)
      end
    end
  end
  return Receptors
end

function SampleReceptor_or_Ligand(State,ParticleType,Species_IDs)

  R = rand(1:sum(State[Species_IDs])) # so we choose the R-th receptor/ligand
  # now we need to find the ID of that particle
  count=0
  RP=0
  for (tx,type) in enumerate(ParticleType)
    if type in Species_IDs
      count+=1
      if count==R
        RP=tx
        break
      end
    end
  end
  #@show count,R
  return RP

end

function sample_ligand_position(pos,rho,h,N)

  # place the recptor at the origin, sample a position from the circle of radius rho. With voxel width h we calculate 
  # the displacement of the voxel that the ligand will lie in.
  wrappedidx = i -> mod(i-1,N) + 1

  r = rho * sqrt(rand())
  theta = 2*pi*rand()
  circle_X = r*cos(theta) 
  circle_Y = r*sin(theta)

  n = maximum([ceil(Int,rho/h-1/2),1])
  gmin = -(n+1/2)*h
  gmax = (n+1/2)*h
  gridpoints = range(gmin,gmax,step=h)
  xdisplacement = 1
  for z=2:length(gridpoints)
      if circle_X > gridpoints[z]
          xdisplacement+=1
      end
  end
  ydisplacement = 1
  for z=2:length(gridpoints)
      if circle_Y > gridpoints[z]
          ydisplacement+=1
      end
  end

  LX = wrappedidx(pos[1] + xdisplacement-(n+1))
  LY = wrappedidx(pos[2] + ydisplacement-(n+1))

  return [LX,LY]
end

function sample_2ligand_binding(PosX,PosY,PlacementRates,RateGridCentre,N)
    wrappedidx = i -> mod(i-1,N) + 1

    PosZ = copy(PosX)
    Jsep = abs(PosX[2]-PosY[2])    # lattice distance between voxels in j index
    Isep = abs(PosX[1]-PosY[1])    # lattice distance between voxels in i index

    # wrapped lengths
    if Jsep > N/2
        Jsep = N - Jsep
    end
    if Isep > N/2
        Isep = N - Isep
    end
    # @show RateGridCentre,Isep,Jsep,PosX,PosY
    Rates = PlacementRates[:,:,RateGridCentre+Isep,RateGridCentre+Jsep]
    SX=findall(x->x>0,Rates)
    Weights=Rates[SX]
    KX=sum(sum(Weights)*rand() .> cumsum(Weights))+1 # index of the voxel that the product is placed into
    if PosY[1]>PosX[1]
        PosZ[1] = wrappedidx(PosX[1] + abs(SX[KX][2]-RateGridCentre))
    else
        PosZ[1] = wrappedidx(PosY[1] + abs(SX[KX][2]-RateGridCentre))
    end
    if PosY[2]>PosX[2]
        PosZ[2] = wrappedidx(PosX[2] + abs(SX[KX][1]-RateGridCentre))
    else
        PosZ[2] = wrappedidx(PosY[2] + abs(SX[KX][1]-RateGridCentre))
    end
    return PosZ
end

function pos2idx(PosX,L,h)
  PosBin = range(0.0,L-h,step=h)

  IX = sum(PosX[1] .> PosBin)
  JX = sum(PosX[2] .> PosBin)

  return [IX,JX]
end

function wrappedpos(PosX,L)

  for i=1:2
      if PosX[i] < 0
          PosX[i] += L
      end
      if PosX[i] > L
          PosX[i] += -L
      end
  end
  return PosX

end

function sample_2ligand_unbinding(PosZ,rho,L,h)
  # Sample the placement of the two molecules X and Y
  # First sample Z within its box
  posZ = h*rand(1,2)
  posZ[1] += h*(PosZ[1]-1)
  posZ[2] += h*(PosZ[2]-1)
  # Next sample X within a ball of radius rho/2
  R = rho/2*sqrt(rand())
  Theta = 2*pi*rand()
  posX = copy(posZ)
  posX[1] += R*cos(Theta)
  posX[2] += R*sin(Theta)
  # Finally align the B molecule through the point of reflection
  posY = 2*posZ-posX
  # Implement periodic boundary conditions
  posX = wrappedpos(posX,L)
  posY = wrappedpos(posY,L)
  # Now find the appropriate indices for the boxes they fall into
  PosX = pos2idx(posX,L,h)
  PosY = pos2idx(posY,L,h)
  return PosX, PosY
end

function implement_reaction!(State,RX,Reactions,SpeciesMap,ParticleType,Pos,cumsum_RN,BoundLigand,Nhbd,Doicutoffd,Gausscutoffd,PlacementRates,DoiRateGridCentre)
  reaction = sum(RX .> cumsum_RN)+1 # tells us which reaction type fired
  # Create two arrays which we can use to optimise update_propensity!()
  reactants = Vector() # stores the types of the particles involved in the reaction
  IDs = Vector() # stores the positions of all particles involved in the reaction

  if ReactionType[reaction]=="Receptor Diffusion"
    # implement a diffusive receptor reaction
    wrappedidx = i -> mod(i-1,N) + 1
    RP = SampleReceptor_or_Ligand(State,ParticleType,ReceptorSpecies)
    # Receptors = findRec_or_Lig(ParticleType,ReceptorSpecies2)
    # RP = rand(Receptors) # choose a random receptor
    dir = rand(1:4) # choose a random direction
    if dir==1
        Pos[RP,1]+=1
        Pos[RP,1] = wrappedidx(Pos[RP,1])
    elseif dir==2
        Pos[RP,1]+=-1
        Pos[RP,1] = wrappedidx(Pos[RP,1])
    elseif dir==3
        Pos[RP,2]+=1
        Pos[RP,2] = wrappedidx(Pos[RP,2])
    else
        Pos[RP,2]+=-1
        Pos[RP,2] = wrappedidx(Pos[RP,2])
    end
    # update the neighbourhood
    update_neighbourhood!(Nhbd,RP,Pos,NTotal,N,Doicutoffd,Gausscutoffd)
    # Store reactants and their IDs
    push!(reactants,ParticleType[RP])
    push!(IDs,RP)
  elseif ReactionType[reaction]=="Ligand Diffusion"
    # implement a diffusive ligand reaction
    wrappedidx = i -> mod(i-1,N) + 1
    RP = SampleReceptor_or_Ligand(State,ParticleType,LigandSpecies)
    # Ligands = findRec_or_Lig(ParticleType,LigandSpecies2)
    # RP = rand(Ligands) # choose a random receptor
    dir = rand(1:4) # choose a random direction
    if dir==1
        Pos[RP,1]+=1
        Pos[RP,1] = wrappedidx(Pos[RP,1])
    elseif dir==2
        Pos[RP,1]+=-1
        Pos[RP,1] = wrappedidx(Pos[RP,1])
    elseif dir==3
        Pos[RP,2]+=1
        Pos[RP,2] = wrappedidx(Pos[RP,2])
    else
        Pos[RP,2]+=-1
        Pos[RP,2] = wrappedidx(Pos[RP,2])
    end
    # update the neighbourhood
    update_neighbourhood!(Nhbd,RP,Pos,NTotal,N,Doicutoffd,Gausscutoffd)
    # Store reactants and their IDs
    push!(reactants,ParticleType[RP])
    push!(IDs,RP)
  else
    reac = Reactions[reaction-2] # gets the data for the reaction
    if ReactionType[reaction] == "First order place none" || ReactionType[reaction] == "First order place one" || ReactionType[reaction] == "First order place two"
      RP = RX - cumsum_RN[reaction-1] # the particular particle that undergoes the reaction
      if ParticleType[RP] != SpeciesMap[reac.substrates[1].op]
        @warn "The species is of type " ParticleType[RP] "it should be of type " SpeciesMap[reac.substrates[1].op]
        @show reaction
      end
      push!(reactants,ParticleType[RP])
      push!(IDs,RP)
      ParticleType[RP] = SpeciesMap[reac.products[1].op]
      push!(reactants,ParticleType[RP])
      push!(IDs,RP)
      if ReactionType[reaction] == "First order place one"
        # now we have to place the ligand that has unbound from the receptor
        LX = BoundLigand[RP]
        BoundLigand[RP]=0
        ParticleType[LX] = SpeciesMap[reac.products[2].op]
        Pos[LX,:] = sample_ligand_position(Pos[RP,:], rho, h, N)
        # update the neighbourhood
        update_neighbourhood!(Nhbd,LX,Pos,NTotal,N,Doicutoffd,Gausscutoffd)
        # store reactants and their IDs
        push!(reactants,ParticleType[LX])
        push!(IDs,LX)

      elseif ReactionType[reaction] == "First order place two"
        # now we have to place the two ligands that unbind from each other
        LX = BoundLigand[RP]
        BoundLigand[RP]=0
        ParticleType[LX] = SpeciesMap[reac.products[2].op]
        Pos[RP,:],Pos[LX,:] = sample_2ligand_unbinding(Pos[RP,:],rho,L,h)
        # update the neighbourhood
        update_neighbourhood!(Nhbd,RP,Pos,NTotal,N,Doicutoffd,Gausscutoffd)
        update_neighbourhood!(Nhbd,LX,Pos,NTotal,N,Doicutoffd,Gausscutoffd)
        # store reactants and their IDs
        push!(reactants,ParticleType[LX])
        push!(IDs,LX)
      end

    elseif ReactionType[reaction] == "Bimolecular Doi"
      Rec = ceil(Int,(RX-cumsum_RN[reaction-1])/NTotal) # the index of the first substrate: a receptor
      Lig = RX-cumsum_RN[reaction-1]-(Rec-1)*NTotal # the index of the second substrate: a ligand
      # store reactants and their IDs
      push!(reactants,ParticleType[Rec])
      push!(IDs,Rec)
      push!(reactants,ParticleType[Lig])
      push!(IDs,Lig)

      ParticleType[Rec] = SpeciesMap[reac.products[1].op] # updates the receptor type
      push!(reactants,ParticleType[Rec])
      push!(IDs,Rec)
      BoundLigand[Rec] = Lig # stores the vector location of the bound ligand (used when unbinding)
      ParticleType[Lig] = -1 # to reflect it is now inactive
    elseif ReactionType[reaction] == "Bimolecular Doi Ligand"
      Lig1 = ceil(Int,(RX-cumsum_RN[reaction-1])/NTotal) # the index of the first substrate: a ligand
      Lig2 = RX-cumsum_RN[reaction-1]-(Lig1-1)*NTotal # the index of the second substrate: a ligand
      # store reactants and their IDs
      push!(reactants,ParticleType[Lig1])
      push!(IDs,Lig1)
      push!(reactants,ParticleType[Lig2])
      push!(IDs,Lig2)

      ParticleType[Lig1] = SpeciesMap[reac.products[1].op] # updates the ligand type
      push!(reactants,ParticleType[Lig1])
      push!(IDs,Lig1)
      #@show Pos[Lig1,:], Pos[Lig2,:]
      Pos[Lig1,:] = sample_2ligand_binding(Pos[Lig1,:],Pos[Lig2,:],PlacementRates,DoiRateGridCentre,N)
      BoundLigand[Lig1] = Lig2 # stores the vector location of the bound ligand (used when unbinding)
      ParticleType[Lig2] = -1 # to reflect it is now inactive

      # update the neighbourhood
      update_neighbourhood!(Nhbd,Lig1,Pos,NTotal,N,Doicutoffd,Gausscutoffd)
    elseif ReactionType[reaction] == "Bimolecular Gauss"
      Rec1 = ceil(Int,(RX-cumsum_RN[reaction-1])/NTotal) # the index of the first substrate: a receptor
      Rec2 = RX-cumsum_RN[reaction-1]-(Rec1-1)*NTotal # the index of the second substrate: a receptor
      # store reactants and their IDs
      push!(reactants,ParticleType[Rec1])
      push!(IDs,Rec1)
      push!(reactants,ParticleType[Rec2])
      push!(IDs,Rec2)

      #ParticleType[Rec1] = SpeciesMap[reac.products[1].op] # updates the receptor type
      ParticleType[Rec2] = SpeciesMap[reac.products[2].op]
      push!(reactants,ParticleType[Rec2])
      push!(IDs,Rec2)

    end

    # Here we update the State vector that keeps track of total species counts
    reac = Reactions[reaction-2]
    sub = reac.substrates
    prod = reac.products
    for s in sub
      State[SpeciesMap[s.op]] += -1
    end
    for p in prod
      State[SpeciesMap[p.op]] += 1
    end

  end
  return reactants,IDs
end

function update_propensity!(reactants,IDs,a0,Propensities,RX,State,SpeciesMap,Reactions,ParameterMap,ReactDependencies,ParticleType,Nhbd,Pos,DoiRateGridCentre,DoiRateMatrix,GaussRateGridCentre,GaussRateMatrix)
  # Here we update the propensities given that reaction RX has just fired

  reaction = sum(RX .> cumsum_RN)+1 # tells us which reaction type fired
  for d in ReactDependencies[reaction]
    
    for (tt,type) in enumerate(reactants)

      # we cycle through the reactants that were affected in the fired reaction and just update the propensity entries for those particles
      # this should result in much faster code
      ID = IDs[tt]

      if ReactionType[d] == "Receptor Diffusion"
        a0 -= Propensities[1]
        Propensities[1] = 4*Dreceptor/h^2 * sum(State[ReceptorSpecies])
        a0 += Propensities[1]
      elseif ReactionType[d] == "Ligand Diffusion"
        a0 -= Propensities[2]
        Propensities[2] = 4*Dligand/h^2 * sum(State[LigandSpecies])
        a0 += Propensities[2]
      elseif ReactionType[d] == "First order place none" || ReactionType[d] == "First order place one" || ReactionType[d] == "First order place two"
        reac = Reactions[d-2] # Reactions is the data structure that does not include diffusive reactions
        # we calculate the propensity rates for a first order reactions

        # turn off the propensity
        a0 -= Propensities[cumsum_RN[d-1]+ID]
        Propensities[cumsum_RN[d-1]+ID] = 0

        if ParticleType[ID] == SpeciesMap[reac.substrates[1].op]
          Propensities[cumsum_RN[d-1]+ID] = ParameterMap[reac.rate.op]
          a0 += Propensities[cumsum_RN[d-1]+ID]
        end
        
      elseif ReactionType[d] == "Bimolecular Doi" || ReactionType[d] == "Bimolecular Doi Ligand"
        reac = Reactions[d-2] # Reactions is the data structure that does not include diffusive reactions
        # we calculate the propensity rates for Bimolecular Doi reactions
        if type == SpeciesMap[reac.substrates[1].op]
          # we search over the ID-th row for neighbours of the 2nd substrate
          for s = 1:NTotal
            a0 -= Propensities[cumsum_RN[d-1]+(ID-1)*NTotal+s]
            Propensities[cumsum_RN[d-1]+(ID-1)*NTotal+s] = 0
          end
          S = findall(x->x==2,Nhbd[ID,:])
          for (sx,s) in enumerate(S)
            # check to see if the molecules are of the right type
            if ParticleType[ID] == SpeciesMap[reac.substrates[1].op] && ParticleType[s]==SpeciesMap[reac.substrates[2].op] 
              xsep = wrappeddist(Pos[ID,1],Pos[s,1],N)
              ysep = wrappeddist(Pos[ID,2],Pos[s,2],N)
              Propensities[cumsum_RN[d-1]+(ID-1)*NTotal+s] = ParameterMap[reac.rate.op] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
              a0 += Propensities[cumsum_RN[d-1]+(ID-1)*NTotal+s]
            end
          end
        elseif type == SpeciesMap[reac.substrates[2].op]
          # we search over ID-th column for neighbours of the 1st substrate
          index = cumsum_RN[d-1]+ID - NTotal
          
          @inbounds for s = 1:NTotal
            index += NTotal
            a0 -= Propensities[index]
            Propensities[index] = 0
          end
          S = findall(x->x==2,Nhbd[:,ID])
          for (sx,s) in enumerate(S)
            # check to see if the molecules are of the right type
            if ParticleType[ID] == SpeciesMap[reac.substrates[2].op] && ParticleType[s]==SpeciesMap[reac.substrates[1].op] 
              xsep = wrappeddist(Pos[ID,1],Pos[s,1],N)
              ysep = wrappeddist(Pos[ID,2],Pos[s,2],N)
              Propensities[cumsum_RN[d-1]+(s-1)*NTotal+ID] = ParameterMap[reac.rate.op] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
              a0 += Propensities[cumsum_RN[d-1]+(s-1)*NTotal+ID]
            end
          end
        end

      elseif ReactionType[d] == "Bimolecular Gauss"
        reac = Reactions[d-2] # Reactions is the data structure that does not include diffusive reactions
        # we calculate the propensity rates for Bimolecular Gauss reactions
        if type == SpeciesMap[reac.substrates[1].op]
          for s = 1:NTotal
            a0 -= Propensities[cumsum_RN[d-1]+(ID-1)*NTotal+s]
            Propensities[cumsum_RN[d-1]+(ID-1)*NTotal+s] = 0
          end
          # we search over the ID-th row for neighbours of the 2nd substrate
          S = findall(x->x==1,Nhbd[ID,:])
          for (sx,s) in enumerate(S)
            # check to see if the molecules are of the right type
            if ParticleType[ID] == SpeciesMap[reac.substrates[1].op] && ParticleType[s]==SpeciesMap[reac.substrates[2].op] 
              xsep = wrappeddist(Pos[ID,1],Pos[s,1],N)
              ysep = wrappeddist(Pos[ID,2],Pos[s,2],N)
              Propensities[cumsum_RN[d-1]+(ID-1)*NTotal+s] = ParameterMap[reac.rate.op] * GaussRateMatrix[GaussRateGridCentre+xsep,GaussRateGridCentre+ysep]
              a0 += Propensities[cumsum_RN[d-1]+(ID-1)*NTotal+s]
            end
          end
        elseif type == SpeciesMap[reac.substrates[2].op]
          for s = 1:NTotal
            a0 -= Propensities[cumsum_RN[d-1]+(s-1)*NTotal+ID]
            Propensities[cumsum_RN[d-1]+(s-1)*NTotal+ID] = 0
          end
          # we search over ID-th column for neighbours of the 1st substrate
          S = findall(x->x==1,Nhbd[:,ID])
          for (sx,s) in enumerate(S)
            # check to see if the molecules are of the right type
            if ParticleType[ID] == SpeciesMap[reac.substrates[2].op] && ParticleType[s]==SpeciesMap[reac.substrates[1].op] 
              xsep = wrappeddist(Pos[ID,1],Pos[s,1],N)
              ysep = wrappeddist(Pos[ID,2],Pos[s,2],N)
              Propensities[cumsum_RN[d-1]+(s-1)*NTotal+ID] = ParameterMap[reac.rate.op] * GaussRateMatrix[GaussRateGridCentre+xsep,GaussRateGridCentre+ysep]
              a0 += Propensities[cumsum_RN[d-1]+(s-1)*NTotal+ID]
            end
          end
        end

      end
    end
  end

  return a0,Propensities

end

function calculate_reaction_numbers(ReactionType,NTotal)

  # First we calculate how long the propensity vector needs to be
  ReactionNumber = zeros(Int64,length(ReactionType))
  for (rx,typ) in enumerate(ReactionType)
    if typ == "Receptor Diffusion" || typ == "Ligand Diffusion"
       ReactionNumber[rx] = 1
    elseif typ == "First order place none" || typ == "First order place one" || typ == "First order place two"
       ReactionNumber[rx] = NTotal 
    elseif typ == "Bimolecular Doi" || typ == "Bimolecular Gauss" || typ == "Bimolecular Doi Ligand"
       ReactionNumber[rx] = NTotal * NTotal
    end
  end
  return cumsum_RN = cumsum(ReactionNumber)
end

function update_neighbourhood!(Nhbd,RP,Pos,NTotal,N,Doicutoffd,Gausscutoffd)
  # update the row and column of the particle RP
  i = RP
  for j in 1:NTotal
    xsep = wrappeddist(Pos[i,1],Pos[j,1],N)
    ysep = wrappeddist(Pos[i,2],Pos[j,2],N)
    if xsep <= Gausscutoffd && ysep <= Gausscutoffd
      Nhbd[i,j]=1
      Nhbd[j,i]=1
      if xsep <= Doicutoffd && ysep <= Doicutoffd
        Nhbd[i,j]=2
        Nhbd[j,i]=2
      end
    else
      Nhbd[i,j]=0
      Nhbd[j,i]=0
    end
  end
  return Nhbd

end

function initialise(State,SpatialState,N)
  Z = State .* SpatialState
  Pos = rand(1:N,sum(Z),2) # sample random initial positions for the particles
  Type = [i for i in 1:length(Z) for j in 1:Z[i]]
  BoundLigand = zeros(Int64,length(Type)) # an array to store the location of the ligands that attach to receptors
  return Pos,Type,BoundLigand
end

function Create_Nhbd_Matrix(Pos,NTotal,Doicutoffd,Gausscutoffd,N)

  # Initialise the neighbourhood matrix 
  Nhbd= zeros(Int,NTotal,NTotal)
  for i=1:NTotal
      for j=1:NTotal
          xsep = wrappeddist(Pos[i,1] ,Pos[j,1],N)
          ysep = wrappeddist(Pos[i,2] ,Pos[j,2],N)
          # if xsep <= Gausscutoffd && ysep <= Gausscutoffd
          #   Nhbd[i,j]=1
          # end
          if xsep <= Gausscutoffd && ysep <= Gausscutoffd
            Nhbd[i,j]=1 # a one is within Gauss reach
            if xsep <= Doicutoffd && ysep <= Doicutoffd
                Nhbd[i,j]=2 # a two is within Doi reach
            end
        end
      end
  end
  return Nhbd
end

function build_dependency_graph(Reactions,ReactionType)
  # Make an array that stores all the reaction dependencies for each reaction
  ReactDependencies = Vector()
  # add receptor diffusive dependencies
  Dependence=Vector()
  for (tx,typ) in enumerate(ReactionType)
    if typ == "Bimolecular Doi" || typ == "Bimolecular Gauss"
      push!(Dependence,tx)
    end
  end
  push!(ReactDependencies,Dependence)
  # add ligand diffusive dependencies
  Dependence=Vector()
  for (tx,typ) in enumerate(ReactionType)
    if typ == "Bimolecular Doi" || typ == "Bimolecular Doi Ligand" || typ == "Bimolecular Gauss"
      push!(Dependence,tx)
    end
  end
  push!(ReactDependencies,Dependence)
  # add non diffusive reaction dependencies
  for (rx,reac) in enumerate(Reactions)
    Ingredients = [reac.substrates; reac.products]
    Ingredients = [ing.op for ing in Ingredients]
    Dependence = Vector()
    # check to see if we add ligand diffusion dependency
    if ReactionType[rx+2] == "Bimolecular Doi" || ReactionType[rx+2] == "Bimolecular Doi Ligand" || ReactionType[rx+2] == "First order place one" || ReactionType[rx+2] == "First order place two"
      push!(Dependence,2) # add ligand diffusion ID which is 2
    end
    for (ix,reac2) in enumerate(Reactions)
      # we search through all reactions in order to build an array of arrays
      Products = reac2.substrates
      Products = [prod.op for prod in Products]
      Overlap = intersect(Ingredients,Products)
      if length(Overlap) > 0
        push!(Dependence,ix+2) # we add two to account for the 2 diffusive reactions
      end
    end
    push!(ReactDependencies,Dependence)
  end
  return ReactDependencies
end


# Build the reaction network that we will use to store all the possible reactions
# NOTE: on unbinding reactions always have the products such that the first is the receptor the 2nd is the ligand (needs to be this way for ligand placement)
rs = @reaction_network begin
    c1, S + X --> SX
    c2, SX --> S + X
    c3, SX --> SXPh
    c4, SXPh --> SX
    c5, SXPh --> S + X
    c6, I + Y --> IYPh
    c7, IYPh --> I + Y
    c8, IYPh + SXPh --> IYPh + SX
    c9, X + Y --> XY
    c10, XY --> X + Y
end c1 c2 c3 c4 c5 c6 c7 c8 c9 c10
# Order the reactions 
Reactions = rs.eqs
# Create an array that details the type of each reaction with diffusive reactions added
ReactionType = ["Receptor Diffusion",
                "Ligand Diffusion",
                "Bimolecular Doi",
                "First order place one",
                "First order place none",
                "First order place none",
                "First order place one",
                "Bimolecular Doi",
                "First order place one",
                "Bimolecular Gauss",
                "Bimolecular Doi Ligand",
                "First order place two"]


# Extract the states and enumerate them from the reaction_network
Species = species(rs) # a list of the species 
SpeciesMap = speciesmap(rs) # a map enumerating the species

# Describes the initial numbers of each species type
InitialState = [10,100,0,0,10,100,0,0] # [S,X,SX,SXPh,I,Y,IYPh,XY] # order of appearance
NLigand = InitialState[2] + InitialState[3] + InitialState[4] + InitialState[6] + InitialState[7] + 2*InitialState[8] # total number of ligand
NReceptors = InitialState[1] + InitialState[3] + InitialState[4] + InitialState[5] + InitialState[7] # total number of receptors
SpatialState = [1,1,1,1,1,1,1,1] # an array assigning vector 1 to species that are resolved spatially and a 0 to those that are well mixed
ReceptorSpecies2 = [1,0,1,1,1,0,1,0] 
ReceptorSpecies = [1,3,4,5,7]
LigandSpecies = [2,6,8]
LigandSpecies2 = [0,1,0,0,0,1,0,1] # this is for ligand unattached to a receptor
NTotal = sum(InitialState[ReceptorSpecies]) + sum(InitialState[LigandSpecies])
cumsum_RN = calculate_reaction_numbers(ReactionType,NTotal)

# Simulation parameters
N = 128 # Number of voxels
#eps = 10.0 # molecular reach in nms of receptor tether
L = 300 # domain length in nms
h = L/N # voxel width in nms
rho = 5.0 # reaction radius in nms of ligand and receptor 
MR = 15.0 # molecular reach

# Reaction Rates
k1 = 100.0 # binding rate of the X ligand
km1 = 0.1 # unbinding rate of the X ligand
k2 = 100.0 # background phosphorylation rate
km2 = 0.1 # background dephosphorylation rate
k3 = 100.0 # binding rate of the Y ligand
km3 = 0.1 # unbinding rate of the Y ligand
k4 = 10.0 / (6.023 * 1e-7) # signal inhibition rate
k5 = 100.0 # colocalisation rate of the X and Y ligand
km5 = 1.0 # unbinding rate of XY

ParamID = rs.ps
Reactions = rs.eqs
p = [k1,km1,k2,km2,km1,k3,km3,k4,k5,km5] # [c1:c10]
Dreceptor = 125.0 # diffusion of receptors
Dligand = 125.0 # diffusion of ligands

ParameterMap = Dict(ParamID[i] => p[i] for i in 1:length(Reactions))

ReactDependencies = build_dependency_graph(Reactions,ReactionType)

function RunSolver()
  t=0
  tf=10.0
  dt = 1e-3
  count = 1
  tSave = range(dt,tf,step=dt)
  PhosphCount = zeros(Int,1,length(tSave))
  State = copy(InitialState) # [S,X,SX,SXPh,I,Y,IYPh] # order of appearance

  # Load the PretabulatedRates

  # Loads Rate Matrices for Bimolecular reactions
  DoiRateMatrix = readdlm("PretabulatedRates/CRDME_DoiRates_7by7_N128.txt")
  DoiRateGridLength = size(DoiRateMatrix)[1]
  DoiRateGridCentre = Int((DoiRateGridLength+1)/2)
  Doicutoffd = DoiRateGridLength - DoiRateGridCentre

  GaussRateMatrix = readdlm("PretabulatedRates/CRDME_Rates_nanom_55by55Grid_L=300_N=128_eps=15.txt")
  GaussRateGridLength = size(GaussRateMatrix)[1]
  GaussRateGridCentre = Int((GaussRateGridLength+1)/2)
  Gausscutoffd = GaussRateGridLength - GaussRateGridCentre

  # Loads the Placement Rates for the placement of the ligand:ligand product Z
  Data = readdlm("PretabulatedRates/PlacementRates_7by7_N128.txt")
  M = size(Data)[1]
  PlacementRates = zeros(M,M,M,M)
  for i=1:M
      for j=1:M
          for u=1:M
              for v=1:M
                  PlacementRates[u,v,i,j] = Data[u,v+(i-1)*M*M+(j-1)*M]
              end
          end
      end
  end

  
  # initialise the positions of the species that are resolved spatially and provides a list of their species types
  Pos,ParticleType,BoundLigand = initialise(State,SpatialState,N)
  Nhbd = Create_Nhbd_Matrix(Pos,NTotal,Doicutoffd,Gausscutoffd,N)
  a0,Propensities = calculate_propensity(State,SpeciesMap,Reactions,ParameterMap,Pos,ParticleType,ReactionType,Nhbd,DoiRateGridCentre,DoiRateMatrix,Doicutoffd,GaussRateGridCentre,GaussRateMatrix,Gausscutoffd)
 #   ReactionCount = zeros(Int,1,54)
  while t < tf
    # update time
    tau = 1/a0*log(1/rand()) 
    t += tau
    #@show t
    # save the output
    while t>tSave[count]
      PhosphCount[count] += State[4]
      #@show t
      count+=1
      if count>length(tSave)
          break
      end
    end

    # sample a reaction
    RX = sample_reaction(a0,Propensities) 
    # ReactionCount[RX]+=1
    # implement the reaction RX
    reactants,IDs=implement_reaction!(State,RX,Reactions,SpeciesMap,ParticleType,Pos,cumsum_RN,BoundLigand,Nhbd,Doicutoffd,Gausscutoffd,PlacementRates,DoiRateGridCentre)
    # update the propensities for reactions that depend on RX
    #Nhbd = Create_Nhbd_Matrix(Pos,NTotal,Doicutoffd,Gausscutoffd,N)
    #a0,Propensities = calculate_propensity(State,SpeciesMap,Reactions,ParameterMap,Pos,ParticleType,ReactionType,Nhbd,DoiRateGridCentre,DoiRateMatrix,Doicutoffd,GaussRateGridCentre,GaussRateMatrix,Gausscutoffd)
    a0,Propensities = update_propensity!(reactants,IDs,a0,Propensities,RX,State,SpeciesMap,Reactions,ParameterMap,ReactDependencies,ParticleType,Nhbd,Pos,DoiRateGridCentre,DoiRateMatrix,GaussRateGridCentre,GaussRateMatrix) 
    #@show a0,t
   
  end
  return PhosphCount
end

# PhosphCount = RunSolver()
# @show PhosphCount

function ManyRealisation(NSim)
  Data = zeros(Float64,1,1000)
#   ReacCount = zeros(Int,1,cumsum_RN[end])
#   MeanEventTime = zeros(Float64,1,NSim)
  for ns = 1:NSim
    PhosphCount = RunSolver()
    Data += PhosphCount
    @show ns
  end
  return Data
end

@profile CD28Count = RunSolver()

     
    


dosave = false

outfname    = "/home/daniel/Documents/MechanismsForEfficientInhibition/Colocalisation_Code/CatalystCode/CaseNeitherTestData.csv"

# Save the data into a csv file
if dosave
  df = DataFrame()
  for i=1
      df[!,Symbol("Realisation ", i)] = CD28Count[i,:]
  end
CSV.write(outfname, df)
end 
