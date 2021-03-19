using Random, Distributions
using CSV, DataFrames, DelimitedFiles

# This code simulates the CRDME model of receptor and ligand dynamics to study receptor inhibition mechanisms
# The reactions incooporated in this script are as follows
# A + X <--> B      or   CD28 + CD80 <--> CD28:CD80
# B <--> C          or   CD28:CD80 <--> CD28:CD80:Ph
# C --> A + X       or   CD28:CD80:ph --> CD28 + CD80
# D + Y <--> E      or   PD1 + PDL1 <--> PD1:PDL1:Ph
# C + E --> B + E   or   CD28:CD80:Ph + PD1:PDL1:Ph --> CD28:CD80 + PD1:PDL1:Ph
# DBWilson 18/03/2021

# Functions to call
function wrappeddist(X1,X2,N)

    if abs(X1-X2) >= ceil(Int,N/2)
        Dist = N - abs(X1-X2)
    else
        Dist = abs(X1-X2)
    end

    Dist 
end

function Calculate_Propensities(Param,PosL,PosR,ReceptorType,LigandType,RV,LV,h,NhbdAX,NhbdDY,NhbdCE,NhbdXY,NStim,NInhib,NStimL,NInhibL,DoiRateMatrix,DoiRateGridCentre,N)
    # The Param vector stores the parameters in the following order
    # Param[1] = Drecept
    # Param[2] = Dligand
    # Param[3] = k1     : A+X --> B
    # Param[4] = km1    : B --> A+X or C --> A+X
    # Param[5] = k2     : B --> C
    # Param[6] = km2    : C --> B
    # Param[7] = k3     : D+Y --> E
    # Param[8] = km3    : E --> D+Y
    # Param[9] = k4     : C+E --> B+E
    # Params[10] = k5   : X+Y --> Z
    # Params[11] = km5  : Z --> X+Y

    Reaction = [1,2,2+NStim*NStimL,2+NStim*NStimL+NInhib*NInhibL,2+NStim*NStimL+NInhib*NInhibL+NStim,2+NStim*NStimL+NInhib*NInhibL+2*NStim,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib+NStim,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib+2*NStim,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib+2*NStim+NStim*NInhib,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib+2*NStim+NStim*NInhib+NStimL*NInhibL,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib+2*NStim+NStim*NInhib+NStimL*NInhibL+(NStimL+NInhibL)]

    Rtypes = length(RV)
    Ltypes = length(LV)

    Propensities = zeros(Reaction[end])
    # Reaction 1 -- receptor diffusion
    Propensities[Reaction[1]] = sum(RV) * 4 * Param[1] / h^2
    #Reaction 2 -- ligand diffusion
    Active = findall(x->x>0,LigandType)
    Propensities[Reaction[2]] = length(Active) * 4 * Param[2] / h^2
    #Reaction 3 -- A+X-->B
    XX = findall(x->x==1, LigandType) # finds the indices of all the X ligands
    for x=1:length(XX)
        XPos = PosL[XX[x],:] # position of the X ligand
        # We search through all the receptors that are in the neighbourhood of the X ligand for A receptors
        AX = findall(x->x==1, NhbdAX[:,XX[x]])
        for a=1:length(AX)
            if ReceptorType[AX[a]]==1
                APos = PosR[AX[a],:] # position of the A receptor
                xsep = wrappeddist(XPos[1],APos[1],N)
                ysep = wrappeddist(XPos[2],APos[2],N)
                Propensities[Reaction[2]+ (XX[x]-1)*NStim + AX[a]] = Param[3] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
            end
        end
    end
    #Reaction 4 -- D+Y-->E
    YX = findall(x->x==2, LigandType) # finds the indices of all the Y ligands
    for y=1:length(YX)
        YPos = PosL[YX[y],:] # position of the Y ligand
        # We search through all the receptors that are in the neighbourhood of the Y ligand for D receptors
        DX = findall(x->x==1, NhbdDY[:,YX[y]-NStimL])
        for d=1:length(DX)
            if ReceptorType[NStim+DX[d]]==4
                DPos = PosR[NStim+DX[d],:] # position of the D receptor
                xsep = wrappeddist(YPos[1],DPos[1],N)
                ysep = wrappeddist(YPos[2],DPos[2],N)
                Propensities[Reaction[3]+ (YX[y]-NStimL-1)*NInhib + DX[d]] = Param[7] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
            end
        end
    end
    #Reaction 5 -- B --> A+X
    for b=1:NStim
        if ReceptorType[b]==2
            Propensities[Reaction[4]+b] = Param[4]
        end
    end
    #Reaction 6 -- C --> A+X
    for c=1:NStim
        if ReceptorType[c]==3
            Propensities[Reaction[5]+c] = Param[4]
        end
    end
    #Reaction 7 -- E --> D+Y
    for e=1:NInhib
        if ReceptorType[e]==5
            Propensities[Reaction[6]+e] = Param[8]
        end
    end
    #Reaction 8 -- B --> C
    for b=1:NStim
        if ReceptorType[b]==2
            Propensities[Reaction[7]+b] = Param[5]
        end
    end
    #Reaction 9 -- C --> B
    for c=1:NStim
        if ReceptorType[c]==3
            Propensities[Reaction[8]+c] = Param[6]
        end
    end
    #Reaction 10 -- C+E --> B+E
    CX = findall(x->x==3, ReceptorType) # finds the indices of all the C molecules
    for c=1:length(CX)
        CPos = PosR[CX[c],:] # position of the C molecule
        # We search through all the E molecules that are in the neighbourhood of the C molecule
        EX = findall(x->x==1, NhbdCE[CX[c],:])
        for e=1:length(EX)
            if ReceptorType[EX[e]+NStim]==5
                EPos = PosR[EX[e]+NStim,:] # position of the E molecule
                xsep = wrappeddist(CPos[1],EPos[1],N)
                ysep = wrappeddist(CPos[2],EPos[2],N)
                Propensities[Reaction[9] + (CX[c]-1)*NInhib + EX[e]] = Param[9] * GaussRateMatrix[GaussRateGridCentre+xsep,GaussRateGridCentre+ysep]
            end
        end
    end
    #Reaction 11 -- X+Y --> Z
    XX = findall(x->x==1, LigandType) # finds the indices of all the X ligands
    for x=1:length(XX)
        XPos = PosL[XX[x],:] # position of the X ligand
        # We search through all the ligands that are in the neighbourhood of the X ligand for Y ligands
        YX = findall(x->x==1, NhbdXY[XX[x],:])
        for y=1:length(YX)
            if LigandType[YX[y]]==2
                YPos = PosL[NStimL+YX[y],:] # position of the Y ligand
                xsep = wrappeddist(XPos[1],YPos[1],N)
                ysep = wrappeddist(XPos[2],YPos[2],N)
                Propensities[Reaction[10]+ (XX[x]-1)*NInhibL + YX[y]] = Param[10] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
            end
        end
    end
    #Reaction 12 -- Z --> X+Y
    for z=1:(NStimL+NInhibL)
        if LigandType[z]==3
            Propensities[Reaction[11]+z] = Param[11]
        end
    end

    # Outputs both the vector of propensities as well as the running total
    return Propensities, sum(Propensities)
end

function Calculate_Reaction(Propensities,a0)
    # We need the sum of Propensities as well as the running total, which we update to avoid recalculating

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
        #@show sum(Propensities), a0, a0*r[2], cumrate
    end
    return RX
end

function Implement_Reaction!(RX,Param,PosR,PosL,NR,NL,NStim,NInhib,NStimL,NInhibL,rho,N,h,ReceptorType,LigandType,BoundLigand,BoundLigand2Ligand,NhbdAX,NhbdDY,NhbdCE,NhbdXY,a0,Propensities,Doicutoffd,DoiRateGridCentre,DoiRateMatrix,Gausscutoffd,GaussRateGridCentre,GaussRateMatrix,PlacementRates,L)

    # Reaction[1] ~ Receptos Diffusion: size 1
    # Reaction[2] ~ Ligand Diffusion: size 1
    # Reaction[3] ~ A+X--> B: size NStim*NStimL
    # Reaction[4] ~ D+Y --> E: size NInhib*NInhibL
    # Reaction[5] ~ B --> A+X: size NStim
    # eventually pass this to the function, rather than recreate it every time
    Reaction = [1,2,2+NStim*NStimL,2+NStim*NStimL+NInhib*NInhibL,2+NStim*NStimL+NInhib*NInhibL+NStim,2+NStim*NStimL+NInhib*NInhibL+2*NStim,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib+NStim,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib+2*NStim,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib+2*NStim+NStim*NInhib,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib+2*NStim+NStim*NInhib+NStimL*NInhibL,2+NStim*NStimL+NInhib*NInhibL+2*NStim+NInhib+2*NStim+NStim*NInhib+NStimL*NInhibL+(NStimL+NInhibL)]

    if RX==Reaction[1]
        # Receptor moves
        wrappedidx = i -> mod(i-1,N) + 1
        RP = rand(1:NR) # choose a random particle
        dir = rand(1:4) # choose a random direction
        if dir==1
            PosR[RP,1]+=1
            PosR[RP,1] = wrappedidx(PosR[RP,1])
        elseif dir==2
            PosR[RP,1]+=-1
            PosR[RP,1] = wrappedidx(PosR[RP,1])
        elseif dir==3
            PosR[RP,2]+=1
            PosR[RP,2] = wrappedidx(PosR[RP,2])
        else
            PosR[RP,2]+=-1
            PosR[RP,2] = wrappedidx(PosR[RP,2])
        end

        # Update propensities and nhbd matrices for Receptor-Ligand interactions
        if ReceptorType[RP] == 1
            # an unbound stimulatory receptor (A) has moved so check the neighbourhood row of Stimulatory ligands (X)
            i = RP # index of stimulatory receptor (A)
            for j=1:NStimL
                a0-=Propensities[ Reaction[2] + (j-1)*NStim + i]; # remove the old propensities from the running total
                # for each stimulatory ligand (X)
                xsep = wrappeddist(PosR[i,1] ,PosL[j,1],N)
                ysep = wrappeddist(PosR[i,2] ,PosL[j,2],N)
                if xsep <= Doicutoffd && ysep <= Doicutoffd 
                    NhbdAX[i,j]=1
                    if LigandType[j] == 1 # checks to see if the ligand is for a stimulatory receptor (X)
                        Propensities[Reaction[2] + (j-1)*NStim + i] = Param[3] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                    else
                        Propensities[Reaction[2] + (j-1)*NStim + i] = 0
                    end
                else
                    NhbdAX[i,j]=0
                    Propensities[Reaction[2] + (j-1)*NStim + i] = 0
                end 
                a0+=Propensities[Reaction[2] + (j-1)*NStim + i]; # store the new propensities
            end

        elseif ReceptorType[RP] == 4
            # an unbound inhibitory receptor (D) has moved so check the neighbourhood row of inhibitory ligands (Y)
            i = RP-NStim # translated index of the inhibitory receptor (D)
            for j=1:NInhibL
                a0-=Propensities[ Reaction[3] + (j-1)*NInhib + i]; # store the old propensities
                # for each inhibitory ligand
                xsep = wrappeddist(PosR[RP,1] ,PosL[j+NStimL,1],N)
                ysep = wrappeddist(PosR[RP,2] ,PosL[j+NStimL,2],N)
                if xsep <= Doicutoffd && ysep <= Doicutoffd 
                    NhbdDY[i,j]=1
                    if LigandType[j+NStimL] == 2 # checks to see if the ligand is for an inhibitory receptor
                        Propensities[Reaction[3] + (j-1)*NInhib + i] = Param[7] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                    else
                        Propensities[Reaction[3] + (j-1)*NInhib + i] = 0
                    end
                else
                    NhbdDY[i,j]=0
                    Propensities[Reaction[3] + (j-1)*NInhib + i] = 0
                end 
                a0+=Propensities[Reaction[3] + (j-1)*NInhib + i]; # store the new propensities
            end

        else
            # a B, C, or E receptor moved
            # We still update neighbourhood matrix incase one of the A,B or D's becomes a C or E later on
            if ReceptorType[RP]==2 || ReceptorType[RP]==3
                # It is a B or C that has moved
                i=RP
                for j=1:NStimL
                    xsep = wrappeddist(PosR[i,1],PosL[j,1],N)
                    ysep = wrappeddist(PosR[i,2],PosL[j,2],N)
                    if xsep <= Doicutoffd && ysep <= Doicutoffd
                        NhbdAX[i,j]=1
                    else
                        NhbdAX[i,j]=0
                    end
                end
            else
                # it is a PD1 (E) that has moved
                i = RP-NStim
                for j=1:NInhibL
                    # for each CD28 molecule 
                    xsep = wrappeddist(PosR[RP,1],PosL[j+NStimL,1],N)
                    ysep = wrappeddist(PosR[RP,2],PosL[j+NStimL,2],N)
                    if xsep <= Doicutoffd && ysep <= Doicutoffd 
                        NhbdDY[i,j]=1
                    else
                        NhbdDY[i,j]=0
                    end 
                end
            end
        end

        # Update propensities and nhbd matrices for Receptor-Receptor interactions
        if ReceptorType[RP] == 3
            # a phosphorylated costimulatory receptor (C) has moved
            i = RP # index of stimulatory receptor (C)
            for j=1:NInhib
                a0-=Propensities[ Reaction[9] + (i-1)*NInhib + j]; # remove the old propensities from the running total
                # for each stimulatory ligand (X)
                xsep = wrappeddist(PosR[i,1] ,PosR[j+NStim,1],N)
                ysep = wrappeddist(PosR[i,2] ,PosR[j+NStim,2],N)
                if xsep <= Gausscutoffd && ysep <= Gausscutoffd 
                    NhbdCE[i,j]=1
                    if ReceptorType[j+NStim] == 5 # checks to see if the coinhibtory receptor is an E
                        Propensities[ Reaction[9] + (i-1)*NInhib + j] = Param[9] * GaussRateMatrix[GaussRateGridCentre+xsep,GaussRateGridCentre+ysep]
                    else
                        Propensities[ Reaction[9] + (i-1)*NInhib + j] = 0
                    end
                else
                    NhbdCE[i,j]=0
                    Propensities[ Reaction[9] + (i-1)*NInhib + j] = 0
                end 
                a0+=Propensities[ Reaction[9] + (i-1)*NInhib + j]; # store the new propensities
            end
        elseif ReceptorType[RP] == 5
            # a phosphorylated coinhibitory receptor (E) has moved
            i = RP-NStim # index of inhibitory receptor (E)
            for j=1:NStim
                a0-=Propensities[ Reaction[9] + (j-1)*NInhib + i]; # remove the old propensities from the running total
                # for each stimulatory ligand (X)
                xsep = wrappeddist(PosR[RP,1] ,PosR[j,1],N)
                ysep = wrappeddist(PosR[RP,2] ,PosR[j,2],N)
                if xsep <= Gausscutoffd && ysep <= Gausscutoffd 
                    NhbdCE[j,i]=1
                    if ReceptorType[j] == 3 # checks to see if the costimulatory receptor is a C
                        Propensities[ Reaction[9] + (j-1)*NInhib + i] = Param[9] * GaussRateMatrix[GaussRateGridCentre+xsep,GaussRateGridCentre+ysep]
                    else
                        Propensities[ Reaction[9] + (j-1)*NInhib + i] = 0
                    end
                else
                    NhbdCE[j,i]=0
                    Propensities[ Reaction[9] + (j-1)*NInhib + i] = 0
                end 
                a0+=Propensities[ Reaction[9] + (j-1)*NInhib + i]; # store the new propensities
            end
        else
            # if either an A,B or D moves still update the neighbourhood matrix
            if ReceptorType[RP]==1 || ReceptorType[RP]==2
                # It is an A or B that has moved
                i=RP
                for j=1:NInhib
                    xsep = wrappeddist(PosR[i,1],PosR[j+NStim,1],N)
                    ysep = wrappeddist(PosR[i,2],PosR[j+NStim,2],N)
                    if xsep <= Gausscutoffd && ysep <= Gausscutoffd
                        NhbdCE[i,j]=1
                    else
                        NhbdCE[i,j]=0
                    end
                end
            else
                # it is a D that has moved
                i = RP-NStim
                for j=1:NStim
                    xsep = wrappeddist(PosR[RP,1],PosR[j,1],N)
                    ysep = wrappeddist(PosR[RP,2],PosR[j,2],N)
                    if xsep <= Gausscutoffd && ysep <= Gausscutoffd 
                        NhbdCE[j,i]=1
                    else
                        NhbdCE[j,i]=0
                    end 
                end
            end
        end


    elseif RX==Reaction[2]
        # Ligand moves
        wrappedidx = i -> mod(i-1,N) + 1
        Active = findall(x->x>0,LigandType)
        RP = rand(Active) # choose a random ligand from those that are active
        dir = rand(1:4) # choose a random direction
        if dir==1
            PosL[RP,1]+=1
            PosL[RP,1] = wrappedidx(PosL[RP,1])
        elseif dir==2
            PosL[RP,1]+=-1
            PosL[RP,1] = wrappedidx(PosL[RP,1])
        elseif dir==3
            PosL[RP,2]+=1
            PosL[RP,2] = wrappedidx(PosL[RP,2])
        else
            PosL[RP,2]+=-1
            PosL[RP,2] = wrappedidx(PosL[RP,2])
        end

        if LigandType[RP] == 1
            # an stimulatory ligand (X) has moved so check the neighbourhood column of Stimulary receptors
            i = RP
            for j=1:NStim
                a0-=Propensities[ Reaction[2] + (i-1)*NStim + j]; # store the old propensities
                # for each stimulatory ligand
                xsep = wrappeddist(PosR[j,1] ,PosL[i,1],N)
                ysep = wrappeddist(PosR[j,2] ,PosL[i,2],N)
                if xsep <= Doicutoffd && ysep <= Doicutoffd 
                    NhbdAX[j,i]=1
                    if ReceptorType[j] == 1 # checks to see if the recepetor is of type A
                        Propensities[Reaction[2] + (i-1)*NStim + j] = Param[3] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                    else
                        Propensities[Reaction[2] + (i-1)*NStim + j] = 0
                    end
                else
                    NhbdAX[j,i]=0
                    Propensities[Reaction[2] + (i-1)*NStim + j] = 0
                end 
                a0+=Propensities[Reaction[2] + (i-1)*NStim + j] # store the new propensities
            end
            # now we check the neighbourhood for the XY reactions
            for j=1:NInhibL
                a0 -= Propensities[Reaction[10] + (i-1)*NInhibL+j]
                # for each Y molecule
                xsep = wrappeddist(PosL[j+NStimL,1],PosL[i,1],N)
                ysep = wrappeddist(PosL[j+NStimL,2],PosL[i,2],N)
                if xsep <= Doicutoffd && ysep <= Doicutoffd 
                    NhbdXY[i,j]=1
                    if LigandType[j+NStimL]==2 # checks to see the Y ligand is active
                        Propensities[Reaction[10] + (i-1)*NInhibL+j] = Param[10] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                    else
                        Propensities[Reaction[10] + (i-1)*NInhibL+j] = 0
                    end
                else
                    NhbdXY[i,j]=0
                    Propensities[Reaction[10] + (i-1)*NInhibL+j] = 0
                end
                a0 +=  Propensities[Reaction[10] + (i-1)*NInhibL+j]
            end

        elseif LigandType[RP] == 2
            # an unbound inhibitory receptor has moved so check the neighbourhood row of inhibitory ligands
            i = RP-NStimL
            for j=1:NInhib
                a0-=Propensities[ Reaction[3] + (i-1)*NInhib + j]; # store the old propensities
                # for each inhibitory ligand
                xsep = wrappeddist(PosR[j+NStim,1] ,PosL[RP,1],N)
                ysep = wrappeddist(PosR[j+NStim,2] ,PosL[RP,2],N)
                if xsep <= Doicutoffd && ysep <= Doicutoffd 
                    NhbdDY[j,i]=1
                    if ReceptorType[j+NStim] == 4 # checks to see if the receptor is of type D
                        Propensities[Reaction[3] + (i-1)*NInhib + j] = Param[7] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                    else
                        Propensities[Reaction[3] + (i-1)*NInhib + j] = 0
                    end
                else
                    NhbdDY[j,i]=0
                    Propensities[Reaction[3] + (i-1)*NInhib + j] = 0
                end 
                a0+=Propensities[Reaction[3] + (i-1)*NInhib + j]; # store the new propensities
            end
            # now we check the neighbourhood for the XY reactions
            for j=1:NStimL
                a0 -= Propensities[Reaction[10] + (j-1)*NInhibL+i]
                # for each Y molecule
                xsep = wrappeddist(PosL[j,1],PosL[i+NStimL,1],N)
                ysep = wrappeddist(PosL[j,2],PosL[i+NStimL,2],N)
                if xsep <= Doicutoffd && ysep <= Doicutoffd 
                    NhbdXY[j,i]=1
                    if LigandType[j]==1 # checks to see the X ligand is active
                        Propensities[Reaction[10] + (j-1)*NInhibL+i] = Param[10] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                    else
                        Propensities[Reaction[10] + (j-1)*NInhibL+i] = 0
                    end
                else
                    NhbdXY[j,i]=0
                    Propensities[Reaction[10] + (j-1)*NInhibL+i] = 0
                end
                a0 +=  Propensities[Reaction[10] + (j-1)*NInhibL+i]
            end

        end

    elseif RX <= Reaction[3]
        # A stimulatory receptor binds to a ligand, i.e. A + X --> B
        XX = ceil(Int,(RX-Reaction[2])/NStim) # gives the index of the X ligand
        AX = RX-Reaction[2]-(XX-1)*NStim # gives the index of the A receptor that becomes a B receptor
        #a0 -= Propensities[1+CX] 
        LigandType[XX]=-1 # this ligand location is now inactive
        BoundLigand[AX]=XX # stores the ligand index for the ligand that is now bound to the receptor 
        Active = findall(x->x>0,LigandType)
        a0-=Propensities[2]
        Propensities[2] = Param[2]/h^2 * length(Active)
        a0+=Propensities[2]
        ReceptorType[AX] = 2
        # turn off the particular A+X--> B reaction 
        a0 -= Propensities[Reaction[2] + (XX-1)*NStim + AX]
        Propensities[Reaction[2] + (XX-1)*NStim + AX] = 0
        # turn on the B --> A+X
        a0-=Propensities[Reaction[4]+AX]
        Propensities[Reaction[4]+AX] = Param[4]
        a0+=Propensities[Reaction[4]+AX]
        # turn on the B --> C
        a0-=Propensities[Reaction[7]+AX]
        Propensities[Reaction[7]+AX] = Param[5]
        a0+=Propensities[Reaction[7]+AX]
        # turn off all other A+X reactions involving the X ligand
        AX2 = findall(x->x==1, NhbdAX[:,XX])
        for a=1:length(AX2)
            if ReceptorType[AX2[a]]==1
                a0 -= Propensities[Reaction[2] + (XX-1)*NStim + AX2[a]]
                Propensities[Reaction[2] + (XX-1)*NStim + AX2[a]] = 0
            end
        end
        # turn off all other A+X reactions involving the A receptor
        XX2 = findall(x->x==1, NhbdAX[AX,:])
        for x=1:length(XX2)
            if LigandType[XX2[x]]==1
                a0 -= Propensities[Reaction[2] + (XX2[x]-1)*NStim + AX]
                Propensities[Reaction[2] + (XX2[x]-1)*NStim + AX] = 0
            end
        end
        # turn off all other X+Y reactions involving the X ligand
        YX= findall(x->x==1, NhbdXY[XX,:])
        for y=1:length(YX)
            if LigandType[YX[y]+NStimL]==2
                a0 -= Propensities[Reaction[10]+(XX-1)*NInhibL + YX[y]]
                Propensities[Reaction[10]+(XX-1)*NInhibL + YX[y]] = 0
            end
        end
        
    elseif RX <= Reaction[4]
        # An inhibitory ligand binds to a receptor, i.e. D + Y --> E
        YX = ceil(Int,(RX-Reaction[3])/NInhib) # gives the index of the Y ligand
        DX = RX-Reaction[3]-(YX-1)*NInhib # gives the index of the D receptor that becomes a E receptor
        #a0 -= Propensities[1+CX] 
        LigandType[YX+NStimL]=-1
        BoundLigand[DX+NStim] = YX
        Active = findall(x->x>0,LigandType)
        a0-=Propensities[2]
        Propensities[2] = Param[2]/h^2 * length(Active)
        a0+=Propensities[2]
        ReceptorType[DX+NStim] = 5
        # turn off particular D+Y-->E reaction
        a0 -= Propensities[Reaction[3] + (YX-1)*NInhib + DX]
        Propensities[Reaction[3] + (YX-1)*NInhib + DX] = 0
        #turn on the E --> D+Y reaction
        a0-=Propensities[Reaction[6]+DX]
        Propensities[Reaction[6]+DX] = Param[8]
        a0+=Propensities[Reaction[6]+DX]
        # turn off all other D+Y reactions involving the Y ligand
        DX2 = findall(x->x==1, NhbdDY[:,YX])
        for d=1:length(DX2)
            if ReceptorType[NStim+DX2[d]]==4
                a0 -= Propensities[Reaction[3] + (YX-1)*NInhib + DX2[d]]
                Propensities[Reaction[3] + (YX-1)*NInhib + DX2[d]] = 0
            end
        end
        # turn off all other A+X reactions involving the A receptor
        YX2 = findall(x->x==1, NhbdDY[DX,:])
        for y=1:length(YX2)
            if LigandType[YX2[y]+NStimL]==2
                a0 -= Propensities[Reaction[3] + (YX2[y]-1)*NInhib + DX]
                Propensities[Reaction[3] + (YX2[y]-1)*NInhib + DX] = 0
            end
        end
        # turn on all C+E-->B+E reactions
        CX = findall(x->x==1, NhbdCE[:,DX])
        for c=1:length(CX)
            if ReceptorType[CX[c]]==3
                a0-=Propensities[Reaction[9] + (CX[c]-1)*NInhib + DX]
                xsep = wrappeddist(PosR[CX[c],1],PosR[DX+NStim,1],N)
                ysep = wrappeddist(PosR[CX[c],2],PosR[DX+NStim,2],N)
                Propensities[Reaction[9] + (CX[c]-1)*NInhib + DX] = Param[9] * GaussRateMatrix[GaussRateGridCentre+xsep,GaussRateGridCentre+ysep]
                a0+=Propensities[Reaction[9] + (CX[c]-1)*NInhib + DX]
            end
        end
        # turn off all other X+Y reactions involving the Y ligand
        XX= findall(x->x==1, NhbdXY[:,YX])
        for x=1:length(XX)
            if LigandType[XX[x]]==1
                a0 -= Propensities[Reaction[10]+(XX[x]-1)*NInhibL + YX]
                Propensities[Reaction[10]+(XX[x]-1)*NInhibL + YX] = 0
            end
        end

    elseif RX <= Reaction[5]
        # a B receptor dissociates from its bound ligand into an A and X
        BX = RX - Reaction[4]
        XX = BoundLigand[BX]
        BoundLigand[BX]=0
        LigandType[XX]=1
        ReceptorType[BX]=1
        a0 -= Propensities[Reaction[4]+BX]
        Propensities[Reaction[4]+BX]=0 # turn off B --> A+X reaction
        Active = findall(x->x>0,LigandType)
        a0-=Propensities[2]
        Propensities[2] = Param[2]/h^2 * length(Active)
        a0+=Propensities[2]
        # Sample a new location for the X ligand
        PosL[XX,:] = sample_ligand_position(PosR[BX,:], rho, h, N)
        # turn on all A+X reactions involving the X ligand
        i = XX
        for j=1:NStim
            a0-=Propensities[ Reaction[2] + (i-1)*NStim + j]; # store the old propensities
            # for each stimulatory ligand
            xsep = wrappeddist(PosR[j,1] ,PosL[i,1],N)
            ysep = wrappeddist(PosR[j,2] ,PosL[i,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdAX[j,i]=1
                if ReceptorType[j] == 1 # checks to see if the recepetor is of type A
                    Propensities[Reaction[2] + (i-1)*NStim + j] = Param[3] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[2] + (i-1)*NStim + j] = 0
                end
            else
                NhbdAX[j,i]=0
                Propensities[Reaction[2] + (i-1)*NStim + j] = 0
            end 
            a0+=Propensities[Reaction[2] + (i-1)*NStim + j] # store the new propensities
        end
        # turn on all A+X reactions involving the A receptor
        i = BX # index of new stimulatory receptor (A)
        for j=1:NStimL
            a0-=Propensities[ Reaction[2] + (j-1)*NStim + i]; # remove the old propensities from the running total
            # for each stimulatory ligand (X)
            xsep = wrappeddist(PosR[i,1] ,PosL[j,1],N)
            ysep = wrappeddist(PosR[i,2] ,PosL[j,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdAX[i,j]=1
                if LigandType[j] == 1 # checks to see if the ligand is for a stimulatory receptor (X)
                    Propensities[Reaction[2] + (j-1)*NStim + i] = Param[3] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[2] + (j-1)*NStim + i] = 0
                end
            else
                NhbdAX[i,j]=0
                Propensities[Reaction[2] + (j-1)*NStim + i] = 0
            end 
            a0+=Propensities[Reaction[2] + (j-1)*NStim + i]; # store the new propensities
        end
        # turn on all X+Y reactions involving the X ligand
        i = XX
        for j=1:NInhibL
            a0-=Propensities[ Reaction[10] + (i-1)*NInhibL + j]; # store the old propensities
            # for each stimulatory ligand
            xsep = wrappeddist(PosL[j+NStimL,1] ,PosL[i,1],N)
            ysep = wrappeddist(PosL[j+NStimL,2] ,PosL[i,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdXY[i,j]=1
                if LigandType[j+NStimL] == 2 # checks to see if the ligand is of type Y
                    Propensities[Reaction[10] + (i-1)*NInhibL + j] = Param[10] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[10] + (i-1)*NInhibL + j] = 0
                end
            else
                NhbdXY[i,j]=0
                Propensities[Reaction[10] + (i-1)*NInhibL + j] = 0
            end 
            a0+=Propensities[Reaction[10] + (i-1)*NInhibL + j] # store the new propensities
        end
        # turn off the B --> C reaction
        a0-=Propensities[Reaction[7]+BX]
        Propensities[Reaction[7]+BX]=0
    elseif RX <= Reaction[6]
        # a C receptor dissociates from its bound ligand into an A and X
        CX = RX - Reaction[5]
        XX = BoundLigand[CX]
        BoundLigand[CX]=0
        LigandType[XX]=1
        ReceptorType[CX]=1
        a0 -= Propensities[Reaction[5]+CX]
        Propensities[Reaction[5]+CX]=0 # turn off B --> A+X reaction
        Active = findall(x->x>0,LigandType)
        a0-=Propensities[2]
        Propensities[2] = Param[2]/h^2 * length(Active)
        a0+=Propensities[2]
        # Sample a new location for the X ligand
        PosL[XX,:] = sample_ligand_position(PosR[CX,:], rho, h, N)
        
        
        
        # remove the C-->B reaction 
        a0 -= Propensities[Reaction[8]+CX]
        Propensities[Reaction[8]+CX] = 0

        # turn off C+E-->B+E reactions
        EX = findall(x->x==1, NhbdCE[CX,:])
        for e=1:length(EX)
            if ReceptorType[EX[e]+NStim]==5
                a0-=Propensities[Reaction[9] + (CX-1)*NInhib + EX[e]]
                Propensities[Reaction[9] + (CX-1)*NInhib + EX[e]] = 0.0
            end
        end
        
        
        # turn on all A+X reactions involving the X ligand
        i = XX
        for j=1:NStim
            a0-=Propensities[ Reaction[2] + (i-1)*NStim + j]; # store the old propensities
            # for each stimulatory ligand
            xsep = wrappeddist(PosR[j,1] ,PosL[i,1],N)
            ysep = wrappeddist(PosR[j,2] ,PosL[i,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdAX[j,i]=1
                if ReceptorType[j] == 1 # checks to see if the recepetor is of type A
                    Propensities[Reaction[2] + (i-1)*NStim + j] = Param[3] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[2] + (i-1)*NStim + j] = 0
                end
            else
                NhbdAX[j,i]=0
                Propensities[Reaction[2] + (i-1)*NStim + j] = 0
            end 
            a0+=Propensities[Reaction[2] + (i-1)*NStim + j] # store the new propensities
        end
        # turn on all A+X reactions involving the A receptor
        i = CX # index of new stimulatory receptor (A)
        for j=1:NStimL
            a0-=Propensities[ Reaction[2] + (j-1)*NStim + i]; # remove the old propensities from the running total
            # for each stimulatory ligand (X)
            xsep = wrappeddist(PosR[i,1] ,PosL[j,1],N)
            ysep = wrappeddist(PosR[i,2] ,PosL[j,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdAX[i,j]=1
                if LigandType[j] == 1 # checks to see if the ligand is for a stimulatory receptor (X)
                    Propensities[Reaction[2] + (j-1)*NStim + i] = Param[3] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[2] + (j-1)*NStim + i] = 0
                end
            else
                NhbdAX[i,j]=0
                Propensities[Reaction[2] + (j-1)*NStim + i] = 0
            end 
            a0+=Propensities[Reaction[2] + (j-1)*NStim + i]; # store the new propensities
        end
        # turn on all X+Y reactions involving the X ligand
        i = XX
        for j=1:NInhibL
            a0-=Propensities[ Reaction[10] + (i-1)*NInhibL + j]; # store the old propensities
            # for each stimulatory ligand
            xsep = wrappeddist(PosL[j+NStimL,1] ,PosL[i,1],N)
            ysep = wrappeddist(PosL[j+NStimL,2] ,PosL[i,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdXY[i,j]=1
                if LigandType[j+NStimL] == 2 # checks to see if the ligand is of type Y
                    Propensities[Reaction[10] + (i-1)*NInhibL + j] = Param[10] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[10] + (i-1)*NInhibL + j] = 0
                end
            else
                NhbdXY[i,j]=0
                Propensities[Reaction[10] + (i-1)*NInhibL + j] = 0
            end 
            a0+=Propensities[Reaction[10] + (i-1)*NInhibL + j] # store the new propensities
        end
    
    elseif RX <= Reaction[7]
        # an E receptor dissociates from its bound ligand into a D and Y
        EX = RX - Reaction[6]
        YX = BoundLigand[EX+NStim]
        BoundLigand[EX+NStim]=0
        LigandType[YX+NStimL]=2
        ReceptorType[EX+NStim]=4
        a0 -= Propensities[Reaction[6]+EX]
        Propensities[Reaction[6]+EX]=0 # turn off E --> D+Y reaction
        Active = findall(x->x>0,LigandType)
        a0-=Propensities[2]
        Propensities[2] = Param[2]/h^2 * length(Active)
        a0+=Propensities[2]
        # Sample a new location for the Y ligand
        PosL[YX+NStimL,:] = sample_ligand_position(PosR[NStim+EX,:], rho, h, N)
        # turn on all D+Y reactions involving the Y ligand
        i = YX
        for j=1:NInhib
            a0-=Propensities[ Reaction[3] + (i-1)*NInhib + j]; # store the old propensities
            # for each stimulatory ligand
            xsep = wrappeddist(PosR[j+NStim,1] ,PosL[i+NStimL,1],N)
            ysep = wrappeddist(PosR[j+NStim,2] ,PosL[i+NStimL,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdDY[j,i]=1
                if ReceptorType[j+NStim] == 4 # checks to see if the recepetor is of type D
                    Propensities[Reaction[3] + (i-1)*NInhib + j] = Param[3] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[3] + (i-1)*NInhib + j] = 0
                end
            else
                NhbdDY[j,i]=0
                Propensities[Reaction[3] + (i-1)*NInhib + j] = 0
            end 
            a0+=Propensities[Reaction[3] + (i-1)*NInhib + j] # store the new propensities
        end
        # turn on all D+Y reactions involving the D receptor
        i = EX # index of new stimulatory receptor (D)
        for j=1:NInhibL
            a0-=Propensities[ Reaction[3] + (j-1)*NInhib + i]; # remove the old propensities from the running total
            # for each stimulatory ligand (X)
            xsep = wrappeddist(PosR[i+NStim,1] ,PosL[j+NStimL,1],N)
            ysep = wrappeddist(PosR[i+NStim,2] ,PosL[j+NStimL,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdDY[i,j]=1
                if LigandType[j+NStimL] == 2 # checks to see if the ligand is for an inhibitory receptor (Y)
                    Propensities[Reaction[3] + (j-1)*NInhib + i] = Param[3] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[3] + (j-1)*NInhib + i] = 0
                end
            else
                NhbdDY[i,j]=0
                Propensities[Reaction[3] + (j-1)*NInhib + i] = 0
            end 
            a0+=Propensities[Reaction[3] + (j-1)*NInhib + i]; # store the new propensities
        end

        # turn off C+E-->B+E reactions
        CX = findall(x->x==1, NhbdCE[:,EX])
        for c=1:length(CX)
            if ReceptorType[CX[c]]==3
                a0-=Propensities[Reaction[9] + (CX[c]-1)*NInhib + EX]
                Propensities[Reaction[9] + (CX[c]-1)*NInhib + EX] = 0.0
            end
        end

        # turn on all X+Y reactions involving the Y ligand
        j = YX
        for i=1:NStimL
            a0-=Propensities[ Reaction[10] + (i-1)*NInhibL + j]; # store the old propensities
            # for each stimulatory ligand
            xsep = wrappeddist(PosL[j+NStimL,1] ,PosL[i,1],N)
            ysep = wrappeddist(PosL[j+NStimL,2] ,PosL[i,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdXY[i,j]=1
                if LigandType[i] == 1 # checks to see if the ligand is of type X
                    Propensities[Reaction[10] + (i-1)*NInhibL + j] = Param[10] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[10] + (i-1)*NInhibL + j] = 0
                end
            else
                NhbdXY[i,j]=0
                Propensities[Reaction[10] + (i-1)*NInhibL + j] = 0
            end 
            a0+=Propensities[Reaction[10] + (i-1)*NInhibL + j] # store the new propensities
        end


    elseif RX <= Reaction[8]
        # a B phosphorylates into a C receptor
        BX = RX - Reaction[7]
        ReceptorType[BX]=3
        # turn off all B --> A+X reactions involving the B receptor amd turn on C --> A+X reactions
        a0-=Propensities[Reaction[4]+BX]
        Propensities[Reaction[4]+BX] = 0.0
        
        a0-=Propensities[Reaction[5]+BX]
        Propensities[Reaction[5]+BX] = Param[4]
        a0+=Propensities[Reaction[5]+BX]
        # turn on C-->B reaction and turn off B-->C reaction
        a0-=Propensities[Reaction[7]+BX]
        Propensities[Reaction[7]+BX]=0.0

        a0-=Propensities[Reaction[8]+BX]
        Propensities[Reaction[8]+BX]=Param[6]
        a0+=Propensities[Reaction[8]+BX]
        # turn on C+E-->B+E reaction
        EX = findall(x->x==1, NhbdCE[BX,:])
        for e=1:length(EX)
            if ReceptorType[EX[e]+NStim]==5
                a0-=Propensities[Reaction[9] + (BX-1)*NInhib + EX[e]]
                xsep = wrappeddist(PosR[BX,1],PosR[EX[e]+NStim,1],N)
                ysep = wrappeddist(PosR[BX,2],PosR[EX[e]+NStim,2],N)
                Propensities[Reaction[9] + (BX-1)*NInhib + EX[e]] = Param[9] * GaussRateMatrix[GaussRateGridCentre+xsep,GaussRateGridCentre+ysep]
                a0+=Propensities[Reaction[9] + (BX-1)*NInhib + EX[e]]
            end
        end
    elseif RX <= Reaction[9]
        # a C dephosphorylates into a B receptor
        CX = RX - Reaction[8]
        ReceptorType[CX]=2
        # turn off all C --> A+X reactions involving the C receptor amd turn on B --> A+X reactions
        a0-=Propensities[Reaction[5]+CX]
        Propensities[Reaction[5]+CX] = 0.0
        
        a0-=Propensities[Reaction[4]+CX]
        Propensities[Reaction[4]+CX] = Param[4]
        a0+=Propensities[Reaction[4]+CX]

        # turn on B-->C reaction and turn off the C-->B reaction
        a0-=Propensities[Reaction[8]+CX]
        Propensities[Reaction[8]+CX]=0.0

        a0-=Propensities[Reaction[7]+CX]
        Propensities[Reaction[7]+CX]=Param[5]
        a0+=Propensities[Reaction[7]+CX]

        # turn off C+E-->B+E reactions
        EX = findall(x->x==1, NhbdCE[CX,:])
        for e=1:length(EX)
            if ReceptorType[EX[e]+NStim]==5
                a0-=Propensities[Reaction[9] + (CX-1)*NInhib + EX[e]]
                Propensities[Reaction[9] + (CX-1)*NInhib + EX[e]] = 0.0
            end
        end



    elseif RX <= Reaction[10]
        # an E receptor inhibits a C receptor: i.e. E+C --> E+B
        CX = ceil(Int,(RX-Reaction[9])/NInhib)
        EX = RX-Reaction[9]-(CX-1)*NInhib
        # change C into a B
        ReceptorType[CX]=2
        # turn off this particular C+E --> B+E reaction
        a0-=Propensities[Reaction[9] + (CX-1)*NInhib+EX]
        Propensities[Reaction[9] + (CX-1)*NInhib+EX]=0.0
        # turn off C --> B
        a0-=Propensities[Reaction[8]+CX]
        Propensities[Reaction[8]+CX]=0.0
        # turn off C --> A+X
        a0-=Propensities[Reaction[5]+CX]
        Propensities[Reaction[5]+CX] = 0.0
        # turn on B --> C
        a0-=Propensities[Reaction[7]+CX]
        Propensities[Reaction[7]+CX]=Param[5]
        a0+=Propensities[Reaction[7]+CX]
        # turn on B --> A+X
        a0-=Propensities[Reaction[4]+CX]
        Propensities[Reaction[4]+CX] = Param[4]
        a0+=Propensities[Reaction[4]+CX]
        # turn off all other C+E --> B+E reactions with any other E receptors
        EX2 = findall(x->x==1, NhbdCE[CX,:])
        for e=1:length(EX2)
            if ReceptorType[EX2[e]+NStim]==5
                a0 -= Propensities[Reaction[9] + (CX-1)*NInhib + EX2[e]]
                Propensities[Reaction[9] + (CX-1)*NInhib + EX2[e]] = 0
            end
        end
    elseif RX <= Reaction[11]
        # Then an X+Y --> Z reaction occurs
        XX = ceil(Int,(RX-Reaction[10])/NInhibL) # gives the index of the Y ligand
        YX = RX-Reaction[10]-(XX-1)*NInhibL # gives the index of the X ligand
        #a0 -= Propensities[1+CX] 
        LigandType[NStimL+YX]=-1 # we turn the Y ligand to be inactive
        LigandType[XX]=3 # We turn the X ligand into the Z ligand
        BoundLigand2Ligand[XX]=YX # stores the Y ligand index in the position of the X ligand index

        # Samples a position for the Z ligand using the midpoint
        PosZ = sampleZ_placement(PosL[XX,:],PosL[YX+NStimL,:],PlacementRates,DoiRateGridCentre,N)
        PosL[XX,:] = PosZ
        
        # updates the ligand diffusion propensity
        Active = findall(x->x>0,LigandType)
        a0-=Propensities[2]
        Propensities[2] = Param[2]/h^2 * length(Active)
        a0+=Propensities[2]
                
        # turn off the particular X+Y--> Z reaction 
        a0 -= Propensities[Reaction[10] + (XX-1)*NInhibL + YX]
        Propensities[Reaction[10] + (XX-1)*NInhibL + YX] = 0
        
        # turn off all other X+Y reactions involving the X ligand
        XX2 = findall(x->x==1, NhbdXY[:,YX])
        for x=1:length(XX2)
            if LigandType[XX2[x]]==1
                a0 -= Propensities[Reaction[10] + (XX2[x]-1)*NInhibL + YX]
                Propensities[Reaction[10] + (XX2[x]-1)*NInhibL + YX] = 0
            end
        end
        # turn off all other X+Y reactions involving the Y ligand
        YX2 = findall(x->x==1, NhbdXY[XX,:])
        for y=1:length(YX2)
            if LigandType[YX2[y]+NStimL]==2
                a0 -= Propensities[Reaction[10] + (XX-1)*NInhibL + YX2[y]]
                Propensities[Reaction[10] + (XX-1)*NInhibL + YX2[y]] = 0
            end
        end
        # turn off all other A+X-->B reactions involving the X ligand
        AX = findall(x->x==1, NhbdAX[:,XX])
        for a=1:length(AX)
            if ReceptorType[AX[a]]==1
                a0 -= Propensities[Reaction[2] + (XX-1)*NStim + AX[a]]
                Propensities[Reaction[2] + (XX-1)*NStim + AX[a]] = 0
            end
        end
        # turn off all other D+Y-->E reactions involving the Y ligand
        DX = findall(x->x==1, NhbdDY[:,YX])
        for d=1:length(DX)
            if ReceptorType[DX[d]+NStim]==4
                a0 -= Propensities[Reaction[3] + (YX-1)*NInhib + DX[d]]
                Propensities[Reaction[3] + (YX-1)*NInhib + DX[d]] = 0
            end
        end
        # turn on the Z-->X+Y reaction
        a0 -= Propensities[Reaction[11]+XX]
        Propensities[Reaction[11]+XX] = Param[11]
        a0 += Propensities[Reaction[11]+XX]
    else
        # Then a Z --> X+Y reaction occurs
        ZX = RX - Reaction[11]
        XX = ZX
        YX = BoundLigand2Ligand[ZX]

        #update ligand types
        LigandType[XX]=1
        LigandType[YX+NStimL]=2

        #update the new positions of the X and Y ligands
        PosX,PosY = sampleXY_placement(PosL[ZX,:],rho,N,L,h)
        PosL[XX,:] = PosX
        PosL[YX+NStimL,:] = PosY

        # update the diffusion of ligands propensity
        Active = findall(x->x>0,LigandType)
        a0-=Propensities[2]
        Propensities[2] = Param[2]/h^2 * length(Active)
        a0+=Propensities[2]

        # turn off the Z --> X+Y reaction
        a0 -= Propensities[Reaction[11]+ZX]
        Propensities[Reaction[11]+ZX] = 0

        # turn on the X+Y --> Z reactions that involve the X ligand
        i = XX
        for j=1:NInhibL
            a0-=Propensities[ Reaction[10] + (i-1)*NInhibL + j]; # store the old propensities
            # for each stimulatory ligand
            xsep = wrappeddist(PosL[j+NStimL,1] ,PosL[i,1],N)
            ysep = wrappeddist(PosL[j+NStimL,2] ,PosL[i,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdXY[i,j]=1
                if LigandType[j+NStimL] == 2 # checks to see if the ligand is of type Y
                    Propensities[Reaction[10] + (i-1)*NInhibL + j] = Param[10] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[10] + (i-1)*NInhibL + j] = 0
                end
            else
                NhbdXY[i,j]=0
                Propensities[Reaction[10] + (i-1)*NInhibL + j] = 0
            end 
            a0+=Propensities[Reaction[10] + (i-1)*NInhibL + j] # store the new propensities
        end

        # turn on the X+Y --> Z reactions that involve the Y ligand
        j = YX
        for i=1:NStimL
            a0-=Propensities[ Reaction[10] + (i-1)*NInhibL + j]; # store the old propensities
            # for each stimulatory ligand
            xsep = wrappeddist(PosL[j+NStimL,1] ,PosL[i,1],N)
            ysep = wrappeddist(PosL[j+NStimL,2] ,PosL[i,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdXY[i,j]=1
                if LigandType[i] == 1 # checks to see if the ligand is of type X
                    Propensities[Reaction[10] + (i-1)*NInhibL + j] = Param[10] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[10] + (i-1)*NInhibL + j] = 0
                end
            else
                NhbdXY[i,j]=0
                Propensities[Reaction[10] + (i-1)*NInhibL + j] = 0
            end 
            a0+=Propensities[Reaction[10] + (i-1)*NInhibL + j] # store the new propensities
        end
        # turn on the A+X --> B reactions that involve the X ligand
        i = XX
        for j=1:NStim
            a0-=Propensities[ Reaction[2] + (i-1)*NStim + j]; # store the old propensities
            # for each stimulatory ligand
            xsep = wrappeddist(PosR[j,1] ,PosL[i,1],N)
            ysep = wrappeddist(PosR[j,2] ,PosL[i,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdAX[j,i]=1
                if ReceptorType[j] == 1 # checks to see if the recepetor is of type A
                    Propensities[Reaction[2] + (i-1)*NStim + j] = Param[3] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[2] + (i-1)*NStim + j] = 0
                end
            else
                NhbdAX[j,i]=0
                Propensities[Reaction[2] + (i-1)*NStim + j] = 0
            end 
            a0+=Propensities[Reaction[2] + (i-1)*NStim + j] # store the new propensities
        end
        # turn on the D+Y --> E reactions that involve the Y ligand
        i = YX
        for j=1:NInhib
            a0-=Propensities[ Reaction[3] + (i-1)*NInhib + j]; # store the old propensities
            # for each stimulatory ligand
            xsep = wrappeddist(PosR[j+NStim,1] ,PosL[i+NStimL,1],N)
            ysep = wrappeddist(PosR[j+NStim,2] ,PosL[i+NStimL,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd 
                NhbdDY[j,i]=1
                if ReceptorType[j+NStim] == 4 # checks to see if the recepetor is of type D
                    Propensities[Reaction[3] + (i-1)*NInhib + j] = Param[3] * DoiRateMatrix[DoiRateGridCentre+xsep,DoiRateGridCentre+ysep]
                else
                    Propensities[Reaction[3] + (i-1)*NInhib + j] = 0
                end
            else
                NhbdDY[j,i]=0
                Propensities[Reaction[3] + (i-1)*NInhib + j] = 0
            end 
            a0+=Propensities[Reaction[3] + (i-1)*NInhib + j] # store the new propensities
        end
    end

    #return a0,PosL,PosR,LigandType,ReceptorType
    return a0
end

function Create_Nhbd_Matrix(PosL,PosR,NStim,NInhib,NStimL,NInhibL,Doicutoffd,Gausscutoffd,N)

    # Initialise the neighbourhood matrix that connects ligands X with receptors A-C
    NhbdAX = zeros(Int,NStim,NStimL)
    for i=1:NStim
        # for each Stimulatory receptor 
        for j=1:NStimL
            # for each Stimulatory ligand 
            xsep = wrappeddist(PosR[i,1] ,PosL[j,1],N)
            ysep = wrappeddist(PosR[i,2] ,PosL[j,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd
                NhbdAX[i,j]=1
            end
        end
    end

    # Initialise the neighbourhood matrix that connects ligands Y with receptors D-E
    NhbdDY = zeros(Int,NInhib,NInhibL)
    for i=1:NInhib
        # for each Inhibitory receptor 
        for j=1:NInhibL
            # for each Inhibitory ligand 
            xsep = wrappeddist(PosR[NStim+i,1] ,PosL[NStimL+j,1],N)
            ysep = wrappeddist(PosR[NStim+i,2] ,PosL[NStimL+j,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd
                NhbdDY[i,j]=1
            end
        end
    end

    # Initialise the neighbourhood matrix that connects receptors C and receptors E
    NhbdCE = zeros(Int,NStim,NInhib)
    for i=1:NStim
        # for each costimulatory receptor
        for j=1:NInhib
            # for each coinhibitory receptor
            xsep = wrappeddist(PosR[i,1] ,PosR[NStim+j,1],N)
            ysep = wrappeddist(PosR[i,2] ,PosR[NStim+j,2],N)
            if xsep <= Gausscutoffd && ysep <= Gausscutoffd
                NhbdCE[i,j]=1
            end
        end
    end

    # Initialise the neighbourhood matrix that connects ligands X and Y
    NhbdXY = zeros(Int,NStimL,NInhibL)
    for i=1:NStimL
        # for each costimulatory ligand
        for j=1:NInhibL
            # for each coinhibitory ligand
            xsep = wrappeddist(PosL[i,1] ,PosL[NStimL+j,1],N)
            ysep = wrappeddist(PosL[i,2] ,PosL[NStimL+j,2],N)
            if xsep <= Doicutoffd && ysep <= Doicutoffd
                NhbdXY[i,j]=1
            end
        end
    end

    return NhbdAX, NhbdDY, NhbdCE, NhbdXY
end

function sample_ligand_position(Pos,rho,h,N)

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

    LX = wrappedidx(Pos[1] + xdisplacement-(n+1))
    LY = wrappedidx(Pos[2] + ydisplacement-(n+1))

    return [LX,LY]
end

function sampleZ_placement(PosX,PosY,PlacementRates,RateGridCentre,N)
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
    #@show RateGridCentre,Isep,Jsep,PosX,PosY
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

function sampleXY_placement(PosZ,rho,N,L,h)
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

function pos2idx(PosX,L,h)
    PosBin = range(0.0,L-h,step=h)

    IX = sum(PosX[1] .> PosBin)
    JX = sum(PosX[2] .> PosBin)

    return [IX,JX]
end

function RunSolver(Nsim,Param,RIC,LIC,N,h,dt,tf,rho,L)

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
    
    tSave = range(dt,tf,step=dt)
    CCount = zeros(Int,1,length(tSave))

    for nr=1:Nsim

        NR = sum(RIC) # total number of receptors
        NL = sum(LIC) # total number of ligand
        NStim = sum(RIC[1:3]) # total number of costimulatory receptors
        NInhib = sum(RIC[4:5]) # total number of coinhibitory receptors
        NStimL = LIC[1] # total number of stimulatory ligand
        NInhibL = LIC[2] # total number of inhibitory ligand

        PosR = rand(1:N,NR,2)
        PosL = rand(1:N,NL,2)
        ReceptorType = zeros(Int,NR)
        count=0
        for i=1:length(RIC)
            ReceptorType[count+1:count+RIC[i]] .= i
            count+=RIC[i]
        end
        LigandType = zeros(Int,NL)
        count=0
        for i=1:length(LIC)
            LigandType[count+1:count+LIC[i]] .= i
            count+=LIC[i]
        end
        # we temporarily assume that there are no B,C, or E receptors initially
        BoundLigand = zeros(Int,NR) # for ligand bound receptors we store the vector location for each of the ligands
        BoundLigand2Ligand = zeros(Int,NStimL) # for the ligand:ligand complex we store the old vector location of the Y ligand
        t=0
        count=1
        NhbdAX, NhbdDY, NhbdCE, NhbdXY = Create_Nhbd_Matrix(PosL,PosR,NStim,NInhib,NStimL,NInhibL,Doicutoffd,Gausscutoffd,N)
        Propensities,a0 = Calculate_Propensities(Param,PosL,PosR,ReceptorType,LigandType,RIC,LIC,h,NhbdAX,NhbdDY,NhbdCE,NhbdXY,NStim,NInhib,NStimL,NInhibL,DoiRateMatrix,DoiRateGridCentre,N)
        ReactionCount = 0*Propensities
        while t < tf

            t += 1/a0 * log(1/rand()) # updates the time  
            Active = findall(x->x>0,LigandType)
            RX = Calculate_Reaction(Propensities,a0)
        
            ReactionCount[RX]+=1
            #a0,PosL,PosR,LigandType,ReceptorType = Implement_Reaction!(RX,PosR,PosL,NR,NL,NStim,NInhib,NStimL,NInhibL,N,ReceptorType,LigandType,NhbdAX,NhbdDY,a0,Propensities,Doicutoffd,DoiRateGridCentre,DoiRateMatrix)
            a0 = Implement_Reaction!(RX,Param,PosR,PosL,NR,NL,NStim,NInhib,NStimL,NInhibL,rho,N,h,ReceptorType,LigandType,BoundLigand,BoundLigand2Ligand,NhbdAX,NhbdDY,NhbdCE,NhbdXY,a0,Propensities,Doicutoffd,DoiRateGridCentre,DoiRateMatrix,Gausscutoffd,GaussRateGridCentre,GaussRateMatrix,PlacementRates,L)

            while t>tSave[count]
                #@show t,sum(Propensities[1:end]),a0
                CX = findall(x->x==3,ReceptorType)
                CCount[count] += length(CX)
                count+=1
                if count>length(tSave)
                    break
                end
                
            end

        end
    end
    return CCount

   
    
end


function VaryPDL1Ligand(k4_star,k5,km5)
    # Simulation parameters
    N = 128 # Number of voxels
    #eps = 10.0 # molecular reach in nms of receptor tether
    L = 300 # domain length in nms
    h = L/N # voxel width in nms
    rho = 5.0 # reaction radius in nms of ligand and receptor 
    MR = 15.0 # molecular reach

    # Model parameters
    Drecept = 1.25 * 1e2 # diffusive coefficient of receptors
    Dligand = 1.25 * 1e2 # diffusive coefficient of ligands used for rate caluclations

    k1 = 100.0 # CD80 ligand binding rate to CD28:      A+X-->B
    km1 = 0.1 # Ligand bound CD28 dissociation rate:   B-->A+X, C-->A+X
    k2 = 100.0 # CD28 phosphorylation rate:               B-->C
    km2 = 0.1 # background dephosphorylation by SHP2:   C-->B
    km3 = 0.1 # PD1 dissociation rate:                 E-->D+Y
    k3 = 100.0 # PD1 ligand binding rate:               D+Y-->E
    k4 = k4_star / (6.023 * 1e-7) # from Sci Adv paper translated to per molecule 6.03 1e-4 / 6.023 1e-7 is approx. 1000: Or Ying's paper 0.1 / 6.023 * 1e-7
    #k5 = 100.0 # CD80 ligand binding rate to PDL1:      X+Y-->Z
    #km5 = 0.1 # unbinding rate of the CD80:PDL1 ligand complex  Z-->X+Y
    Param = [Drecept,Dligand,k1,km1,k2,km2,k3,km3,k4,k5,km5]

    NCD28 = 10 # number of CD28 receptor molecules. 
    NCD80 = 100 # number of CD80 ligand molecules
    NPD1 = 100 # number of PD1 receptor molecules
    Nsim = 1
  
    TSpan = [1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,100.0,100.0] 
    DTSpan = [1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,1e-2,1e-3,1e-3] 
    PDL1Span = [1, 3, 5, 10, 22, 47, 100, 216, 465, 1000]

    for z = 1:length(PDL1Span)

        tf = TSpan[z]
        dt = DTSpan[z]
        NPDL1 = PDL1Span[z]
  
        # Initial particle numbers: A,B,C,D,E,X,Y
        RIC = [NCD28,0,0,NPD1,0] 
        LIC = [NCD80,NPDL1]

        if RIC[2] > 0 || RIC[3] > 0  || RIC[5] > 0
            @warn "Starting with B, C, or E molecules is not supported, as the vector location for the bound ligand needs to be defined"
        end

        save = true
        @time  CCount = RunSolver(1,Param,RIC,LIC,N,h,dt,tf,rho,L)
        #@show size(CCount)
        # Save the data into a csv file
        if save
            outfname = "SpatialLigandData/SpatialLigandModel-CaseNeither-k4=$(k4_star)-k5=$(k5)-km5=$(km5)-CD80Number=$(NCD80)-PDL1Number=$(NPDL1)-MolecularReach=$(MR)-N=$(N)-PD1=$(NPD1)-CD28=$(NCD28)-tf=$(tf)-dt=$(dt)-nsims=$(Nsim).csv"
            df = DataFrame()
            for i=1:Nsim
                df[!,Symbol("Realisation ", i)] = CCount[i,:]
            end
            CSV.write(outfname, df)
        end 
    end

end

k5 = 100.0
km5 = 0.001
k4_star = 0.0

VaryPDL1Ligand(k4_star,k5,km5)