using Random, Distributions
using CSV, DataFrames, DelimitedFiles

# This code simulates the enhances direct Gillespie method for the 3 species model in the comp biol paper
# Note my indexing does require that we store the PD1 molcules and then the CD28 molecules. They cannot become mixed

N = 64 # number of voxels
# Nsim = 100
# tf = 0.5
# dt = 1e-3

MR = 15 # moecular reach CHANGE THE INPUTS RATE

NCD28 = 10 # NUmber of CD28 receptors
NCD80 = 100 # number of CD80 ligand molecules. NB this just rescales on rate, no ligand interactions are modelled here.
NPD1 = 1 # number of PD1 receptors

NsimSpan = [1,1,1,1,1,1,1,1,1,1]
#I0Span = [1, 5, 22, 100, 465, 2155, 10000, 46416, 215444, 1000000]
PDL1Span = [1, 3, 5, 10, 22, 47, 100, 216, 465, 1000]
#I0Span = [1000000];
TSpan = [1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,100.0,100.0]*100
DTSpan = [1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,1e-3,1e-3]*100

for z = 1:length(PDL1Span)

    Nsim = NsimSpan[z]
    tf = TSpan[z]
    dt = DTSpan[z]
    NPDL1 = PDL1Span[z]

# initial particle numbers: CD28 --> A      CD28-CD80 --> B      CD28-CD80-Ph --> C     PD1 --> D     PD1-PDL2--Ph --> E
# Start with all A and E such that we dont have more ligand bound receptors than there are ligands.
NA = NCD28 # initial number of A molecules
NB = 0 # initial number of B molecules
NC = 0 # initial number of C molecules
ND = NPD1 # initial number of D molecules
NE = 0 # initial number of E molecules

k5_star = 100.0 # phosphorylation rate

    # saving parameters
    dosave      = true
    outfname    = "/home/daniel/Documents/MechanismsForEfficientInhibition/ModellingLigandInteractions/DoseResponseCurves/AllPathways/Ligand_Dimerisation/data/WellMixedLigandModel-CasePD1-k5=$(k5_star)-CD80Number=$(NCD80)-PDL1Number=$(PDL1Span[z])-MolecularReach=$(MR)-N=$(N)-PD1=$(NPD1)-CD28=$(NCD28)-tf=$(tf)-dt=$(dt)-nsims=$(Nsim).csv"

    function wrappeddist(X1,X2,N)

        if abs(X1-X2) >= ceil(Int,N/2)
            Dist = N - abs(X1-X2)
        else
            Dist = abs(X1-X2)
        end

        Dist 
    end

    function RunSolver(Nsim,tf,dt,NA,NB,NC,ND,NE,NCD80,N,wrappeddist,k5_star)
        
        NS = ceil(Int,tf/dt)+1

        CCopyNumber = zeros(NS,Nsim)
        Tspan = 0:dt:(NS-1)*dt

        # model parameters
        eps = 15 # molecular reach in nms
        L = 300 # domain length in nms
        h = L/N # voxel width in nms
        rho = 5 # reaction radius in nm of ligand and receptor, used to calculate the ligand binding rate with the smoluchowski rate

        Drecept = 1.25 * 1e2 # diffusive coefficient of receptors
        Dligand = 1.25 * 1e2 # diffusive coefficient of ligands used for rate caluclations

        k1_star = 100.0 #  CD80 ligand binding rate
        k3_star = 100.0 # PDL1 ligand binding rate
        #k5_star = 100.0 # CD80 to PDL1 association rate
        k1 = k1_star * (25*pi)/(300^2) 
        km1 = 0.1 # CD28 dissociation rate
        k2 = 100.0 # CD28 phosphorylation rate
        km2 = 0.1 # background dephosphorylation by SHP2
        km3 = 0.1 # PD1 dissociation rate
        k3 = k3_star * (25*pi)/(300^2)  # PD1 ligand binding rate
        k4 = 10.0 / (6.023 * 1e-7) # from Sci Adv paper translated to per molecule 6.03 1e-4 / 6.023 1e-7 is approx. 1000: Or Ying's paper 0.1 / 6.023 * 1e-7
        k5 = k5_star * (25*pi)/(300^2) # binding rate of ligand dimerisation
        km5 = 1.0 # unbinding rate of dimerised ligand

        #k4 = 1.492691965991955*1e7 # catalytic efficiency of phosphorylated PD1 and phosphorylated CD28: chosen such that beta=10
        # conserved numbers of inhibitory and activating receptors
        NT = NA+NB+NC+ND+NE # total number of receptors
        NPD1 = ND+NE # total number of inhibitory receptors
        NCD28 = NA+NB+NC # total number of activating receptors

        #Sigma = 7.443673151765994 * 1e-7
        #k4 = Betaspan[y]*k2/(Sigma*NPD1) # catalytic efficiency of phosphorylated PD1 and phosphorylated CD28: chosen to fix beta

        # Rates for bimolecular gaussian reactions
        RateMatrix = readdlm("CRDME_Rates_nanom_31by31Grid_L=300_N=64_eps=15.txt")
        RateGridLength = size(RateMatrix)[1]
        RateGridCentre = Int((RateGridLength+1)/2)
        cutoffd = RateGridLength - RateGridCentre

        for nsim=1:Nsim

            @show z,nsim

            CCopyNumber[1,nsim] = NC # add the initial number of C molcules -- i.e. phosphorylated CD28-CD80 complex

            # Initialise position and type vectors
            Pos = rand(1:N,NT,2)
            ParticleType = zeros(Int,NT,1)
            ParticleType[1:NA].=1
            ParticleType[NA+1:NA+NB].=2
            ParticleType[NA+NB+1:NA+NB+NC].=3
            ParticleType[NA+NB+NC+1:NA+NB+NC+ND].=4
            ParticleType[NA+NB+NC+ND+1:NT].=5

            # Initialise a ligand type. This tells you what type of ligand is bound to B,C and E receptors
            # 1 -- CD80, 2 -- PDL1, 3 -- CD80:PDL1
            BoundLigandType = zeros(Int,NT,1)
            BoundLigandType[NA+1:NA+NB+NC].=1 # initally all B and C receptors are attached to just a CD80
            BoundLigandType[NA+NB+NC+ND+1:NT].=2 # initally all E receptores are attached to just a PDL1 

            current_CD80 = NCD80 - (NB+NC)
            current_PDL1 = NPDL1 - NE
            current_CD80_PDL1 = 0 # initially no bound ligand complexes

            # Initialise the neighbourhood matrix that connects receptors that are close enough to react
            Neighbourhood = zeros(Int,NPD1,NCD28)

            for i=1:NCD28
                # for each CD28 molecule 
                for j=1:NPD1
                    # for each PD1 molecule 
                    xsep = wrappeddist(Pos[i,1] ,Pos[NCD28+j,1],N)
                    ysep = wrappeddist(Pos[i,2] ,Pos[NCD28+j,2],N)
                    if xsep <= cutoffd && ysep <= cutoffd
                        Neighbourhood[j,i]=1
                    end
                end
            end
            
            # Store a vector of reaction types
            NReac = 3 + NT+ NPD1*NCD28 # 1 diffusive, NT first order reactions, NCD28*NPD1 bimolecular reactions
            ReactType = zeros(Int,NReac,1)
            ReactType[1] = 1 # these are diffusive reactions 
            ReactType[2] = 2 # these are the ligand binding reactions
            ReactType[3] = 3 # these are the ligand unbinding reactions
            ReactType[4:3+NT] .= 4 # these are all the other first order reactions involving  receptors and/or ligands
            ReactType[3+NT+1:end] .= 5 # these are second order inhibition reactions

            # Initialise reactions vector
            Propensities = zeros(Float64,NReac)
            # Add the diffusive propensities for the all molecules
            Propensities[1] = NT*4*Drecept/h^2 
            # Add the ligand dimerisation reaction
            Propensities[2] = k5*current_CD80*current_PDL1
            # Add the ligand de-dimerisation reaction
            Propensities[3] = km5*current_CD80_PDL1
            # Add the first order propensities
            Propensities[4:3+NT] = k1*(current_CD80)*(ParticleType .== 1) + (k2+km1)*(ParticleType .== 2) + (km2+km1)*(ParticleType .== 3) + k3*(current_PDL1+current_CD80_PDL1)*(ParticleType .== 4) + km3*(ParticleType .== 5)
            # Add the bimolecular dephosphorylation propensities
            CX = findall(x->x==3, ParticleType) # finds the indices of all the C molecules
            for c=1:length(CX)
                CPos = Pos[CX[c],:] # position of the C molecule
                # We search through all the E molecules that are in the neighbourhood of the C molecule
                EX = findall(x->x==1, Neighbourhood[:,CX[c]])
                for e=1:length(EX)
                    EPos = Pos[EX[e]+NCD28,:] # position of the E molecule
                    xsep = wrappeddist(CPos[1],EPos[1],N)
                    ysep = wrappeddist(CPos[2],EPos[2],N)
                    Propensities[3+NT+ (CX[c]-1)*NPD1 + EX[e]] = k4 * RateMatrix[RateGridCentre+xsep,RateGridCentre+ysep]
                end
            end


            # now we intialise the algorithm
            wrappedidx = i -> mod(i-1,N) + 1
            t = 0.0 # current time
            a0 = sum(Propensities) # total propensity rate
            count = 2
            RX=0
            ss=1.0
            while t < tf
                OldRX = RX
                OldPropensities = Propensities
                r = rand(2) # two random numbers to select the next reaction time and the next reaction
                t += 1/a0 * log(1/r[1]) # updates the time 
                #@show current_CD80, current_PDL1, current_CD80_PDL1, Propensities[2], Propensities[3]
                #sleep(0.01)
                while t>Tspan[count] && count<NS
                    # we save the total number of C molecules
                    CCopyNumber[count,nsim] = length(findall(x->x==3, ParticleType))
                    count+=1
                end

                if t < tf # no need to simulate the reaction if it occurs after tf
                    if t > ss
                    #@show t
                    ss+=1.0
                    end
                    # calculate the cumulative sum as you go
                    reacted=0
                    cumrate=0
                    RX=1
                    while reacted==0
                        if RX==4+NT+NPD1*NCD28
                            @show a0*r[2], a0, r[2], sum(Propensities)
                        end
                        cumrate += Propensities[RX]
                        if a0*r[2] < cumrate
                            reacted=1
                        else
                        RX+=1
                        end
                        #@show sum(Propensities), a0, a0*r[2], cumrate
                    end

                    #RX = sum(a0*r[2] .> cumsum(Propensities)) + 1 # selects the reaction that occurs

                    if ReactType[RX] == 1
                        # implement a diffusive reaction and update the necessary propensities
                        RP = rand(1:NT) # choose a random particle
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
                        # Only update bimolecular propensities after jump if the jump was a molecule of type C or E
                        if ParticleType[RP] == 5
                            # a PD1 receptor has moved so check the neighbourhood row
                            j = RP-NCD28
                            for i=1:NCD28
                                a0-=Propensities[3+NT + (i-1)*NPD1 + j]; # store the old propensities
                                # for each CD28 molecule 
                                xsep = wrappeddist(Pos[i,1] ,Pos[RP,1],N)
                                ysep = wrappeddist(Pos[i,2] ,Pos[RP,2],N)
                                if xsep <= cutoffd && ysep <= cutoffd 
                                    Neighbourhood[j,i]=1
                                    if ParticleType[i] == 3 # checks to see if the CD28 is of type C
                                        Propensities[3+NT + (i-1)*NPD1 + j] = k4 * RateMatrix[RateGridCentre+xsep,RateGridCentre+ysep]
                                    else
                                        Propensities[3+NT + (i-1)*NPD1 + j] = 0
                                    end
                                else
                                    Neighbourhood[j,i]=0
                                    Propensities[3+NT + (i-1)*NPD1 + j] = 0
                                end 
                                a0+=Propensities[3+NT + (i-1)*NPD1 + j]; # store the new propensities
                            end
                        elseif ParticleType[RP]==3
                            # a CD28 receptor has moved so check the neighbourhood column
                            i = RP
                            for j=1:NPD1
                                a0-=Propensities[3+NT + (i-1)*NPD1 + j]; # store the old propensities
                                # for each PD1 molecule 
                                xsep = wrappeddist(Pos[i,1],Pos[j+NCD28,1],N)
                                ysep = wrappeddist(Pos[i,2],Pos[j+NCD28,2],N)
                                if xsep <= cutoffd && ysep <= cutoffd
                                    Neighbourhood[j,i]=1
                                    if ParticleType[j+NCD28] == 5
                                        Propensities[3+NT + (i-1)*NPD1 + j] = k4 * RateMatrix[RateGridCentre+xsep,RateGridCentre+ysep]
                                    else
                                        Propensities[3+NT + (i-1)*NPD1 + j] = 0
                                    end
                                else
                                    Neighbourhood[j,i]=0
                                    Propensities[3+NT + (i-1)*NPD1 + j] = 0
                                end
                                a0+=Propensities[3+NT + (i-1)*NPD1 + j]; # store the new propensities
                            end
                        else
                            # We still update neighbourhood matrix incase one of the A,B or D's becomes a C or E later on
                            if ParticleType[RP]<3
                                # It is a CD28 (A or B) that has moved
                                i=RP
                                for j=1:NPD1
                                    xsep = wrappeddist(Pos[i,1],Pos[j+NCD28,1],N)
                                    ysep = wrappeddist(Pos[i,2],Pos[j+NCD28,2],N)
                                    if xsep <= cutoffd && ysep <= cutoffd
                                        Neighbourhood[j,i]=1
                                    else
                                        Neighbourhood[j,i]=0
                                    end
                                end
                            else
                                # it is a PD1 (D) that has moved
                                j = RP-NCD28
                                for i=1:NCD28
                                    # for each CD28 molecule 
                                    xsep = wrappeddist(Pos[i,1],Pos[RP,1],N)
                                    ysep = wrappeddist(Pos[i,2],Pos[RP,2],N)
                                    if xsep <= cutoffd && ysep <= cutoffd 
                                        Neighbourhood[j,i]=1
                                    else
                                        Neighbourhood[j,i]=0
                                    end 
                                end
                            end
                        end
                    elseif ReactType[RX] == 2
                        if current_PDL1==0
                            @warn "This cant have happened"
                        end
                        # a CD80 and PDL1 ligand dimerise
                        current_CD80+=-1
                        current_PDL1+=-1
                        current_CD80_PDL1+=1

                        # update the relevant propensities
                        a0 -= Propensities[2]
                        a0 -= Propensities[3]
                        Propensities[2] = k5*current_CD80*current_PDL1
                        Propensities[3] = km5*current_CD80_PDL1
                        a0 += Propensities[2]
                        a0 += Propensities[3]
                        for k=1:NCD28
                            if ParticleType[k]==1
                                a0 -= Propensities[k+3]
                                Propensities[k+3] = k1*(current_CD80)
                                a0 += Propensities[k+3]
                            end
                        end
                        for k=1:NPD1
                            if ParticleType[k+NCD28]==4
                                a0 -= Propensities[3+NCD28+k]
                                Propensities[3+NCD28+k] = k3 * (current_PDL1+current_CD80_PDL1)
                                a0 += Propensities[3+NCD28+k]
                            end
                        end


                    elseif ReactType[RX] == 3
                        current_CD80+=1
                        current_PDL1+=1
                        current_CD80_PDL1+=-1

                        # update the relevant propensities
                        a0 -= Propensities[2]
                        a0 -= Propensities[3]
                        Propensities[2] = k5*current_CD80*current_PDL1
                        Propensities[3] = km5*current_CD80_PDL1
                        a0 += Propensities[2]
                        a0 += Propensities[3]
                        for k=1:NCD28
                            if ParticleType[k]==1
                                a0 -= Propensities[k+3]
                                Propensities[k+3] = k1*(current_CD80)
                                a0 += Propensities[k+3]
                            end
                        end
                        for k=1:NPD1
                            if ParticleType[k+NCD28]==4
                                a0 -= Propensities[3+NCD28+k]
                                Propensities[3+NCD28+k] = k3 * (current_PDL1+current_CD80_PDL1)
                                a0 += Propensities[3+NCD28+k]
                            end
                        end

                    elseif ReactType[RX] == 4       
                        # now we have to check which first order reaction occurs
                        RP = RX-3 # this gives th index of the reactant particle, should lie between 1 and NT
                        if ParticleType[RP]==1
                            # then CD28 binds to a ligand, i.e. A --> B
                            ParticleType[RP]=2
                            current_CD80-=1
                               
                            a0 -= Propensities[RX]
                            Propensities[RX] = (k2+km1)
                            a0 += Propensities[RX]
                            # A ligand has been used so update all other relevant propensities
                            for k=1:NCD28
                                if ParticleType[k]==1
                                    a0 -= Propensities[k+3]
                                    Propensities[k+3] = k1*(current_CD80)
                                    a0 += Propensities[k+3]
                                end
                            end
                            
                            a0 -= Propensities[2]
                            Propensities[2] = k5*current_CD80*current_PDL1
                            a0 += Propensities[2]
                            
                        elseif ParticleType[RP]==2
                            # then bound CD28 either phosphorylates or dissociates, i.e. B --> C or B --> A
                            runif = rand(1)
                            if runif[1] < km1/(k2+km1)
                                # then bound CD28 dissociates
                                ParticleType[RP]=1
                                current_CD80 += 1
                                
                                # A ligand has been used so update all other relevant propensities
                                for k=1:NCD28
                                    if ParticleType[k]==1
                                        a0 -= Propensities[k+3]
                                        Propensities[k+3] = k1*(current_CD80)
                                        a0 += Propensities[k+3]
                                    end
                                end
                                
                                a0 -= Propensities[2]
                                Propensities[2] = k5*current_CD80*current_PDL1
                                a0 += Propensities[2]
                                
                            else
                                # then bound CD28 phosphorylates
                                ParticleType[RP]=3
                                a0 -= Propensities[RX]
                                Propensities[RX]=km1+km2
                                a0 += Propensities[RX]
                                
                                # turn on the necessary bimolecular reactions
                                CPos = Pos[RP,:]
                                EX = findall(x->x==1,Neighbourhood[:,RP])
                                for e=1:length(EX)
                                    if ParticleType[EX[e]+NCD28]==5
                                        EPos = Pos[EX[e]+NCD28,:] # position of the E molecule
                                        xsep = wrappeddist(CPos[1],EPos[1],N)
                                        ysep = wrappeddist(CPos[2],EPos[2],N)
                                        Propensities[3+NT+(RP-1)*NPD1 + EX[e]] = k4 * RateMatrix[RateGridCentre+xsep,RateGridCentre+ysep]
                                        a0 += Propensities[3+NT+(RP-1)*NPD1 + EX[e]]
                                    end
                                end
                                
                            end
                        elseif ParticleType[RP]==3
                            runif = rand(1)
                            if runif[1] < km2/(km2+km1)
                                # then phosphorylated CD28 dephoshorylates, i.e. C --> B
                                ParticleType[RP]=2
                                a0 -= Propensities[RX]
                                Propensities[RX]=(k2+km1)
                                a0 += Propensities[RX]
                                # turn off the necessary bimolecular reactions
                                EX = findall(x->x==1,Neighbourhood[:,RP])
                                for e=1:length(EX)
                                    a0 -= Propensities[3+NT+(RP-1)*NPD1 + EX[e]]
                                    Propensities[3+NT+(RP-1)*NPD1 + EX[e]] = 0
                                end

                            else
                                # then phosphorylated CD28 dissociates, i.e. C --> A
                                ParticleType[RP]=1
                                current_CD80 += 1
                                
                                # A ligand has been involved so update all other relevant propensities
                                for k=1:NCD28
                                    if ParticleType[k]==1
                                        a0 -= Propensities[k+3]
                                        Propensities[k+3] = k1 * (current_CD80)
                                        a0 += Propensities[k+3]
                                    end
                                end
                                
                                a0 -= Propensities[2]
                                Propensities[2] = k5*current_CD80*current_PDL1
                                a0 += Propensities[2]
                                # turn off the necessary bimolecular reactions
                                EX = findall(x->x==1,Neighbourhood[:,RP])
                                for e=1:length(EX)
                                    a0 -= Propensities[3+NT+(RP-1)*NPD1 + EX[e]]
                                    Propensities[3+NT+(RP-1)*NPD1 + EX[e]] = 0
                                end
                            end
                        
                        elseif ParticleType[RP]==4
                            # then PD1 binds to ligand, i.e. D--> E
                            ParticleType[RP]=5
                            if current_PDL1==0 && current_CD80_PDL1==0
                                @warn "This reaction shouldn't have happened"
                                sleep(5.0)
                            end

                            if rand() < current_PDL1/(current_PDL1+current_CD80_PDL1)
                                #then D binds to a PDL1
                                current_PDL1 -= 1
                                BoundLigandType[RP]=2
                            else
                                # D binds to a CD80:PDL1 complex
                                current_CD80_PDL1 -= 1
                                BoundLigandType[RP]=3
                            end
                            a0 -= Propensities[RX]
                            Propensities[RX]=km3
                            a0 += Propensities[RX]

                            a0 -= Propensities[2]
                            Propensities[2] = k5*current_CD80*current_PDL1
                            a0 += Propensities[2]
                            # A ligand has been involved so update all other relevant propensities
                            for k=1:NPD1
                                if ParticleType[k+NCD28]==4
                                    a0 -= Propensities[3+NCD28+k]
                                    Propensities[3+NCD28+k] = k3 * (current_PDL1+current_CD80_PDL1)
                                    a0 += Propensities[3+NCD28+k]
                                end
                            end
                            
                            # turn on the necessary bimolecular reactions
                            EPos = Pos[RP,:] # position of the E molcule
                            CX = findall(x->x==1,Neighbourhood[RP-NCD28,:])
                            for c=1:length(CX)
                                if ParticleType[CX[c]]==3
                                    CPos = Pos[CX[c],:] # position of the C molecule
                                    xsep = wrappeddist(CPos[1],EPos[1],N)
                                    ysep = wrappeddist(CPos[2],EPos[2],N)
                                    Propensities[3+NT+(CX[c]-1)*NPD1 + RP-NCD28] = k4 * RateMatrix[RateGridCentre+xsep,RateGridCentre+ysep]
                                    a0 += Propensities[3+NT+(CX[c]-1)*NPD1 + RP-NCD28]
                                end
                            end
                            
                        elseif ParticleType[RP]==5
                            # then phosphorylated PD1 dissociates, i.e. E--> D
                            ParticleType[RP]=4
                            if BoundLigandType[RP]==2
                                current_PDL1 += 1
                            elseif BoundLigandType[RP]==3
                                current_CD80_PDL1 += 1
                            end
                            BoundLigandType[RP]=0
                            # A ligand has been involved so update all other relevant propensities
                            for k=1:NPD1
                                if ParticleType[k+NCD28]==4
                                    a0 -= Propensities[3+NCD28+k]
                                    Propensities[3+NCD28+k] = k3 * (current_PDL1+current_CD80_PDL1)
                                    a0 += Propensities[3+NCD28+k]
                                end
                            end
                            
                            a0 -= Propensities[2]
                            Propensities[2] = k5*current_CD80*current_PDL1
                            a0 += Propensities[2]
                            # turn off the necessary bimolecular reactions
                            CX = findall(x->x==1,Neighbourhood[RP-NCD28,:])
                            for c=1:length(CX)
                                a0 -= Propensities[3+NT+(CX[c]-1)*NPD1 + RP-NCD28]
                                Propensities[3+NT+(CX[c]-1)*NPD1 + RP-NCD28] = 0
                            end
                            
                        end

                    elseif ReactType[RX] == 5
                        # implement a dephosphorylation reaction and update the necessary propensities
                        CX = ceil(Int,(RX-3-NT)/NPD1) # gives the indice of the C molecule that becomes a B
                        a0 -= Propensities[3+CX] 
                        ParticleType[CX] = 2
                        Propensities[3+CX] = (k2+km1)
                        a0 += Propensities[3+CX]
                        EX = findall(x->x==1, Neighbourhood[:,CX])
                        for e=1:length(EX)
                            a0 -= Propensities[3+NT + (CX-1)*NPD1 + EX[e]]
                            Propensities[3+NT + (CX-1)*NPD1 + EX[e]] = 0
                        end
                        
                    end
                    
                 #@show (a0 - sum(Propensities))/(sum(Propensities))

                end

                
                
            end
            CCopyNumber[NS,nsim] = length(findall(x->x==3, ParticleType))
            #NC = length(findall(x->x==3, ParticleType))
            #NB = length(findall(x->x==2, ParticleType))
            #NE = length(findall(x->x==5, ParticleType))
            #@show current_CD80+current_CD80_PDL1+NB+NC, current_PDL1+current_CD80_PDL1+NE

        end
        CCopyNumber
    end


    @time CCopyNumber = RunSolver(Nsim,tf,dt,NA,NB,NC,ND,NE,NCD80,N,wrappeddist,k5_star)

    # Save the data into a csv file
    if dosave
            df = DataFrame()
            for i=1:Nsim
                df[!,Symbol("Realisation ", i)] = CCopyNumber[:,i]
            end
        CSV.write(outfname, df)
    end 
    
end
