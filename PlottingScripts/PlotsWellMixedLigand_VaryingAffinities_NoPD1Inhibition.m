%% This code plots the CD28-Ph fraction as a function of ligand 

close all
clear all
clc


% NPD1Span = 100;

Str1a = '../WellMixedLigandData/WellMixedLigandModel-CaseNeither-k4=0.0-k5=100.0-km5=';
Str3 = '-CD80Number=';
Str5 = '-PDL1Number=';
Str7a = '-MolecularReach=15-N=64-PD1=';
Str9 = '-CD28=10-tf=';
Str11 = '-dt=';
Str13 = '-nsims=1.csv';

Km5Span = [0.001,0.01,0.1,1.0,10.0,100.0];
for km5s=1:length(Km5Span)
    km5 = Km5Span(km5s);

    NPD1 = 100;
    NCD28=10;
    NCD80=100;
    if km5 >= 1.0
        Str2 = num2str(km5,'%.1f');
    else
        Str2 = num2str(km5);
    end

    Str4 = num2str(NCD80);

    %     I0Span = [1, 5, 22, 100, 465, 2155, 10000, 46416, 215444, 1000000];
    PDL1Span = [1, 3, 5, 10, 22, 47, 100, 216, 465, 1000];
    %I0Span = [1, 4, 13, 47];
    %     I0Span(10)=[];
    if NPD1==1
        TSpan = [1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,100.0,100.0]*100;
        DTSpan = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.1,0.1];
    elseif NPD1==1000
        TSpan = [1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,100.0,100.0];
        DTSpan = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.1,0.1]/100;
    else
        TSpan = [1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,100.0,100.0];
        DTSpan = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.1,0.1]/100;
    end


    CD28Ph = zeros(1,length(PDL1Span));
    for z=1:10

        Str6 = int2str(PDL1Span(z));

        Str8 = int2str(NPD1);
        Str10 = num2str(TSpan(z),'%.1f');
        if NPD1==1
            Str12 = num2str(DTSpan(z),'%.1f');
        else
            Str12 = num2str(DTSpan(z));
        end

        T = readtable(strcat(Str1a,Str2,Str3,Str4,Str5,Str6,Str7a,Str8,Str9,Str10,Str11,Str12,Str13));
        T = table2array(T);
%         plot(T)
%         pause
        %         figure(3)
        %         plot(T)10
        %         pause
        %         hold on
        CD28Ph(z) = mean(T(10000:end));

    end

    semilogx(PDL1Span,CD28Ph/10,'*-','MarkerSize',10)
    hold on
end

%% Run the case with no ligand binding to each other

Str1a = '../WellMixedLigandData/WellMixedLigandModel-CaseNeither-k4=0.0-k5=0.0-km5=';
Str3 = '-CD80Number=';
Str5 = '-PDL1Number=';
Str7a = '-MolecularReach=15-N=64-PD1=';
Str9 = '-CD28=10-tf=';
Str11 = '-dt=';
Str13 = '-nsims=1.csv';

Km5Span = [1.0];
for km5s=1:length(Km5Span)
    km5 = Km5Span(km5s);

    NPD1 = 100;
    NCD28=10;
    NCD80=100;
    if km5 >= 1.0
        Str2 = num2str(km5,'%.1f');
    else
        Str2 = num2str(km5);
    end

    Str4 = num2str(NCD80);

    %     I0Span = [1, 5, 22, 100, 465, 2155, 10000, 46416, 215444, 1000000];
    PDL1Span = [1, 3, 5, 10, 22, 47, 100, 216, 465, 1000];
    %I0Span = [1, 4, 13, 47];
    %     I0Span(10)=[];
    if NPD1==1
        TSpan = [1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,100.0,100.0]*100;
        DTSpan = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.1,0.1];
    elseif NPD1==1000
        TSpan = [1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,100.0,100.0];
        DTSpan = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.1,0.1]/100;
    else
        TSpan = [1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,1000.0,100.0,100.0];
        DTSpan = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.1,0.1]/100;
    end


    CD28Ph = zeros(1,length(PDL1Span));
    for z=1:10

        Str6 = int2str(PDL1Span(z));

        Str8 = int2str(NPD1);
        Str10 = num2str(TSpan(z),'%.1f');
        if NPD1==1
            Str12 = num2str(DTSpan(z),'%.1f');
        else
            Str12 = num2str(DTSpan(z));
        end

        T = readtable(strcat(Str1a,Str2,Str3,Str4,Str5,Str6,Str7a,Str8,Str9,Str10,Str11,Str12,Str13));
        T = table2array(T);
        %         figure(3)
        %         plot(T)10
        %         pause
        %         hold on
        CD28Ph(z) = mean(T(10000:end));

    end

    semilogx(PDL1Span,CD28Ph/10,'dk-','MarkerSize',10)
    hold on
end

%% Add figure text
legend('km5 = 0.001','km5 = 0.01','km5 = 0.1','km5 = 1.0','km5 = 10.0','km5 = 100.0','No dimerisation')
xlabel('Number of PDL1 ligand')
ylabel('fraction of phosphorylated CD28')
title('Well mixed ligand dimerisation: case neither')
