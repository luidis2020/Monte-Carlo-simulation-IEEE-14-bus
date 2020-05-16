%--------------------------------------------------------------------------
%--------------------- Programa de Fluxo de Potência ----------------------
%--------------------------------------------------------------------------
%clc;
%clear all;
    
%pdf_load_ieee14;
MVAbase=100;
itmax=10;
control=0;
tol=0.0001;
marcafluxo = 1;
plotgraf=1;
%deg=2;
%  deg = menu (' Deseja dar um degrau na Carga?    ',...
%            '                  Sim                   ',...
%            '                  Não                  ',...
%            '                Cancela                ');
%        
% for deg = 3
%    % k=0;
%         break
% end
       
% if deg == 1
%     
%     deg1 = menu ('Faça sua Escola!  ',...
%                  'P = 1.05 e Q = 1.05  ',...
%                  '      P = 1.00 e Q = 1.05    ',...
%                  '      P = 1.05 e Q = 1.00    ',...
%                  '      P = 0.95 e Q = 0.95    ',...
%                  '      P = 1.05 e Q = 0.95    ',... 
%                  '      P = 0.95 e Q = 1.05    ',...
%                  '            Cancela          '); 
%     if deg1 == 1
%        dgp = 1.05;
%        dgq = 1.05;
%     end
%     if deg1 == 2
%        dgp = 1.00;
%        dgq = 1.05;
%     end
% 
%     if deg1 == 3
%        dgp = 1.05;
%        dgq = 1.00;
%     end
%     if deg1 == 4
%      dgp = 0.95;
%      dgq = 0.95;
%     end
%     
%     if deg1 == 5
%       dgp = 1.05;
%       dgq = 0.95;
%     end
% 
%     if deg1 == 6
%        dgp = 0.95;
%        dgq = 1.05;
%     end
% 
%     if deg1 == 7
%        k = 0;
% %        break
%     end
% 
% end

% if deg == 2
%      dgp=degrau_gp(MCS_iter);
%      %dgq = ;
%      deg1 = 0;
% end

%--------------------------------------------------------------------------
% --------------- Inialização das Variáveis do Programa -------------------
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% -------------------- Leitura dos Dados do Sistema -----------------------
%--------------------------------------------------------------------------
k=1;
%    k = menu('              Escolha o Sistema Teste!         ',...
%             '              Two area symmetric system        ',...
%             '                      Cancela                  ');

if k==1
   %Leitura dos Dados dos arquivos armazenados em branch - Kundur
   %[numb,nomebus,narea,nzl,tipobus,V,T,pc,qc,pg,qg,base,vdes,qmax,qmin,gshk,bshk,rem]=...
   %textread('pai_barras_ipfc.txt','%f%7c%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
   %[de,para,arearamo,ls,circ,tiporamo,r,x,bshl,nl1,nl2,nl3,cont,lado,tap,fi,...
   %tapmin,tapmax,ssize,tmax,tmin]=...
   %textread('pai_ramos_ipfc.txt','%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f');
  
   eval('case_IEEE14_general')
   %eval('two_area_symmetric_original')
   
   sistema = 1;
   marcador = 1;   
end
for k=2
    break
end

%MCS_IEEE14;

if  marcador == 1 

    for l = 1:1
    numb  = barras_ipfc(:,1); 
    nomebus = barras_ipfc(:,2); 
    narea =  barras_ipfc(:,3); 
    nzl    = barras_ipfc(:,4); 
    tipobus  = barras_ipfc(:,5); 
    V   = barras_ipfc(:,6); 
    T   = barras_ipfc(:,7); 
    
%   pc   = demanda_PG(:,1); 
%    qc   = demanda_PG(:,2); 
%     
     pc   = barras_ipfc(:,8); 
     qc   = barras_ipfc(:,9); 
    pg  = barras_ipfc(:,10); 
    qg  = barras_ipfc(:,11); 
    base = barras_ipfc(:,12); 
    vdes = barras_ipfc(:,13); 
    qmax = barras_ipfc(:,14); 
    qmin = barras_ipfc(:,15); 
    gshk = barras_ipfc(:,16); 
    bshk = barras_ipfc(:,17); 
    rem = barras_ipfc(:,18); 
    end

    for l = 1:1
    de  = ramos_ipfc(:,1); 
    para = ramos_ipfc(:,2); 
    arearamo =  ramos_ipfc(:,3); 
    ls    = ramos_ipfc(:,4); 
    circ  = ramos_ipfc(:,5); 
    tiporamo  = ramos_ipfc(:,6); 
    r   = ramos_ipfc(:,7); 
    x   = ramos_ipfc(:,8); 
    bshl   = ramos_ipfc(:,9); 
    nl1  = ramos_ipfc(:,10); 
    nl2  = ramos_ipfc(:,11); 
    nl3 = ramos_ipfc(:,12); 
    cont = ramos_ipfc(:,13); 
    lado = ramos_ipfc(:,14); 
    tap = ramos_ipfc(:,15); 
    fi = ramos_ipfc(:,16);
    tapmin = ramos_ipfc(:,17);
    tapmax = ramos_ipfc(:,18);
    ssize = ramos_ipfc(:,19);
    tmax = ramos_ipfc(:,20);
    tmin = ramos_ipfc(:,21);
    end

    
    %--------------------------------------------------------------------------
    %------------------- Determinando o Tamanho do Sistema --------------------
    %--------------------------------------------------------------------------

    nb = length(numb); % Número de Barras
    nr = length(de); % Número de Ramos
    limQ(nb,1:2) = 0; %Verificação Limite de Reativo
    radgraus = (180/pi); %Radianos para graus
    grausrad = (pi/180); %Graus para radianos

    %--------------------------------------------------------------------------
    %---------------- Determinando Quantidades do Sistema ---------------------
    %--------------------------------------------------------------------------

    %     Tipo 3 --> slack        Tipo 2 --> pv         Tipo 0 --> pq

    nslack = length(find(tipobus == 3)); % números barras slack
    npv = length(find(tipobus == 2)); % números barras pv
    npq = length(find(tipobus == 0)); % números barras pq
    ntransf = length(find(tap ~=1)); % números de transformadores
    ncargas = length(find(pc ~=0));% números de cargas do sistema
    nshunts = length(find(bshk ~=0)); % números de shunts de barra

    %---------------------- Fim da Quantidade do Sistema ----------------------

    %--------------------------------------------------------------------------
    %-------------------- Convertendo Grandezas em p.u. -----------------------
    %--------------------------------------------------------------------------

    pg = pg/MVAbase; %potência ativa gerada dada em MW
    qg = qg/MVAbase; %potência reativa gerada dada em MVar
    pc = pc/MVAbase; %potência ativa consumida dada em MW
    qc = qc/MVAbase; %potência reativa consumida dada em MVar
    qmax = qmax/MVAbase; %potência reativa máxima dada em MVar
    qmin = qmin/MVAbase; %potência reativa mínima dada em MVar

    %------------------ Fim da Conversão para p.u. ----------------------------

    %--------------------------------------------------------------------------
    %------------------ Algumas Conversões Necessárias ------------------------
    %--------------------------------------------------------------------------

    T = T*grausrad; % Passando de graus para radianos
    fi = fi*grausrad;
    bshl = bshl/2; % Dividindo Shunt da Linha por 2

    %--------------------------------------------------------------------------
    %------------ Cálculo da Impedância e Admitância do Sistema ---------------
    %--------------------------------------------------------------------------

    tipo = tipobus;
    Vesp = V;
    z = r + j*x;
    ykm = 1./z;
    gkm=real(ykm);
    bkm=imag(ykm);

    %----------- Fim do Cálculo da Impedância e Admitância do Sistema ---------

    %--------------------------------------------------------------------------
    %----- Cálculo da Potência Ativa e Reativa Especificada em cada Barra -----
    %--------------------------------------------------------------------------

    Pesp = pg - pc;
    Qesp = qg - qc;

    %-------- Fim do Cálculo da potência ativa e reativa especificada ---------

    %--------------------------------------------------------------------------
    %------------------ Construção da Matriz Admitância -----------------------
    %--------------------------------------------------------------------------

    Y = zeros(nb,nb);

    G = zeros(nb,nb); % Matriz condutância 
    B = zeros(nb,nb); % Matriz susceptância  

    for k=1:nb 
        B(k,k) = bshk(k); %Elementos das diagonais matriz B
    end

    for k=1:nb 
        G(k,k) = gshk(k); %Elementos das diagonais matriz G
    end

    for l=1:nr
        k = de(l);
        m = para(l);

        gkmf = r(l)/(r(l)^2+x(l)^2);
        bkmf = -x(l)/(r(l)^2+x(l)^2);

        % Elementos das diagonais
        G(k,k) = G(k,k)+tap(l)^2*gkmf;
        B(k,k) = B(k,k)+tap(l)^2*bkmf + bshl(l); 

        G(m,m) = G(m,m)+gkmf;
        B(m,m) = B(m,m)+bkmf+bshl(l);

        % Elementos fora das diagonais.
        G(k,m) = G(k,m)-tap(l)*gkmf; 
        B(k,m) = B(k,m)-tap(l)*bkmf;

        G(m,k) = G(m,k)-tap(l)*gkmf;
        B(m,k) = B(m,k)-tap(l)*bkmf;   
    end

    Y = G + j*B;

    %----------------Fim da construção da matriz admitância ---------------

    %----------------------------------------------------------------------
    %- Construção da Matriz Admitância sem shunt da barra nem shunt da
    %linha
    %----------------------------------------------------------------------

    YY = zeros(nb,nb);

    GG = zeros(nb,nb); % Matriz condutância 
    BB = zeros(nb,nb); % Matriz susceptância  

    for k=1:nb 
        BB(k,k) = 0; %Elementos das diagonais matriz B
    end

    for k=1:nb 
        GG(k,k) = 0; %Elementos das diagonais matriz G
    end

    for l=1:nr
        k = de(l);
        m = para(l);

        ggkmf = r(l)/(r(l)^2+x(l)^2);
        bbkmf = -x(l)/(r(l)^2+x(l)^2);

        % Elementos das diagonais
        GG(k,k) = GG(k,k)+tap(l)^2*ggkmf;
        BB(k,k) = BB(k,k)+tap(l)^2*bbkmf; 

        GG(m,m) = GG(m,m)+ggkmf;
        BB(m,m) = BB(m,m)+bbkmf;

        % Elementos fora das diagonais.
        GG(k,m) = GG(k,m)-tap(l)*ggkmf; 
        BB(k,m) = BB(k,m)-tap(l)*bbkmf;

        GG(m,k) = GG(m,k)-tap(l)*ggkmf;
        BB(m,k) = BB(m,k)-tap(l)*bbkmf;   
    end

    YY = GG + j*BB;

    %----------------Fim da construção da matriz admitância -------------------

    %--------------------------------------------------------------------------
    %------- Separando a parte real e a parte imaginaria da Matriz [Y] --------
    %--------------------------------------------------------------------------

    MAP_DELTA_V=[-BB G;-GG -B];
    G = real(Y);
    B = imag(Y);

    for k=1:nb 
        if tipo(k)<=2, % Barra PV ou PQ 
            T(k) = 0;
            if tipo(k)==0, % Barra PQ 
                V(k) = 1.0;
            end
        end
    end

    %--------------------------------------------------------------------------

    %--------------------------------------------------------------------------
    %--------------- Configurações Gerais do processo Iterativo ---------------
    %--------------------------------------------------------------------------

    flag = 0;  %bandeira do contador
    it = 0;    %contador de iterações

    %--------------------------------------------------------------------------

    while (it < itmax && flag==0)

    Pcal = zeros(nb,1);
    Qcal = zeros(nb,1);

    % Cálculo da Potência Líquida

    for k = 1:nb    
     for m = 1:nb

       Tkm = T(k) - T(m);
       Pcal(k)= Pcal(k) + V(k)*V(m)*(G(k,m)*cos(Tkm) + B(k,m)*sin(Tkm));
       Qcal(k)= Qcal(k) + V(k)*V(m)*(G(k,m)*sin(Tkm) - B(k,m)*cos(Tkm));

     end
    end

    %--------------------------------------------------------------------------
    %----- Checando Possíveis Violações de Potência Reativa dos Geradores -----
    %--------------------------------------------------------------------------

    if (control ==1)
    Qg = Qcal + qc;

    if it>0
    for k = 1:nb
        if tipo(k) == 2 %% Verificaçao do limite Q se for barra PV.
           if Qg(k) > qmax(k)   %% Verificando limite superior.
              Qg(k) = qmax(k);
              Qesp(k) = Qg(k) - qc(k);
              tipo(k) = 3;
              limQ(k,2) = 1;
              fprintf('Na Iter %u a Barra %u Violou o Limite Superior\n\n',it,numb(k,:));
              fprintf('----------------------------------------------\n\n');
            elseif Qg(k) < qmin(k)   %% Verificando limite inferior.
                   Qg(k) = qmin(k);
                   tipo(k) = 3;
                   Qesp(k) = Qg(k) - qc(k);
                   limQ(k,1) = 1;
                  fprintf('Na Iter %u a Barra %u Violou o Limite Inferior\n\n',it,numb(k,:));
                  fprintf('----------------------------------------------\n\n');
            end
        end
    end
    end
    end

    %--------------------------- Fim da Verificação ---------------------------


    %--------------------------------------------------------------------------
    %------------------ Cálculo dos Mismatches de Potência --------------------
    %--------------------------------------------------------------------------

    DP = zeros(1,nb);
    DQ = zeros(1,nb);
    DP = Pesp - Pcal;
    DQ = Qesp - Qcal;


    %-------------- Fim do Cálculo dos Mismatches de Potência -----------------

    %--------------------------------------------------------------------------
    % Zerando Mismatches de Potência ativa e reativa nas Barras Slack e reativa
    %---------------------------- nas barras PV -------------------------------
    %--------------------------------------------------------------------------

    for ii = 1: nb
        if (tipo(ii) == 3)
            DP(ii) = 0;
            DQ(ii) = 0;
        elseif (tipo(ii) == 2)
            DQ(ii) = 0;
        end
    end

    %--------------------------------------------------------------------------

    %--------------------------------------------------------------------------         
    %---------------------------- Parte Dinâmica ------------------------------
    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------
    %--------- Criando o Vetor Mismatches de Potência ativa e reativa ---------
    %--------------------------------------------------------------------------

    DPQ = [DP;DQ]; 

    %--------------------------------------------------------------------------
    %----------------- Checando a convergência do Vetor Mismatches ------------
    %--------------------------------------------------------------------------

    for ii = 1:((nb*2))      %Verificar (Deixar variáveis do IPFC fora)
        if (abs(DPQ) < tol)
           flag = 1;
        end
    end

    %--------------------------------------------------------------------------

    %--------------------------------------------------------------------------
    %-------------------- Construção da Matriz Jacobiana ----------------------
    %--------------------------------------------------------------------------

    %--------------------------------------------------------------------------
    %----------------- Criando as submatrizes H N M e L nulas -----------------       
    %--------------------------------------------------------------------------
    %%

    delp = max(abs(DP));
    delq = max(abs(DQ));

    %% Rotina para obter a matriz Jacobiana ------------------------------------
    %%

        H = zeros(nb,nb);
        M = H;
        L = H;
        N = H;
        HESSIAN_H=zeros(nb,nb);
        HESSIAN_M=HESSIAN_H;
        HESSIAN_L=HESSIAN_H;
        HESSIAN_N=HESSIAN_H;

        for k=1:nb 
            H(k,k) = -Qcal(k)-V(k)*V(k)*B(k,k); %Elementos da diagonal 
            HESSIAN_H(k,k)=-Pcal(k);
            N(k,k) = (Pcal(k)+V(k)*V(k)*G(k,k))/V(k);
            HESSIAN_N(k,k)=2*G(k,k);
            M(k,k) = Pcal(k)-V(k)*V(k)*G(k,k);
            HESSIAN_M(k,k)=-Qcal(k);
            L(k,k) = (Qcal(k)-V(k)*V(k)*B(k,k))/V(k); 
            HESSIAN_L(k,k)=-2*B(k,k);
        end

        for l=1:nr 
            k  = de(l);
            m  = para(l);
            ab = T(k) - T(m); % Ângulo entre a barra K e a barra M. 

            gkmf = tap(l)*r(l)/(r(l)^2+x(l)^2); % Condutâncias especificadas entre as barras
            bkmf = -tap(l)*x(l)/(r(l)^2+x(l)^2);% Susceptâncias especificadas entre as barras

            % Expressões das derivadas 
            H(k,m) = H(k,m)+V(k)*V(m)*(-gkmf*sin(ab)+bkmf*cos(ab));
            HESSIAN_H(k,m)=HESSIAN_H(k,m)-V(k)*V(m)*(gkmf*cos(ab)+bkmf*sin(ab));
            H(m,k) = H(m,k)-V(k)*V(m)*(-gkmf*sin(ab)-bkmf*cos(ab));
            HESSIAN_H(m,k)=HESSIAN_H(m,k)-V(k)*V(m)*(gkmf*cos(ab)-bkmf*sin(ab));
            %%%%%%%%%%%%%%
            N(k,m) = N(k,m)+ V(k)*(-gkmf*cos(ab)-bkmf*sin(ab));
            HESSIAN_N(k,m)=HESSIAN_N(k,m)+(-gkmf*cos(ab)-bkmf*sin(ab));
            N(m,k) = N(m,k)+ V(m)*(-gkmf*cos(ab)+bkmf*sin(ab));
            HESSIAN_N(m,k)=HESSIAN_N(m,k)+(-gkmf*cos(ab)+bkmf*sin(ab));
            %%%%%%%%%%%%%%
            M(k,m) = M(k,m)-V(k)*V(m)*(-gkmf*cos(ab)-bkmf*sin(ab));
            HESSIAN_M(k,m)=HESSIAN_M(k,m)-V(k)*V(m)*(+gkmf*sin(ab)-bkmf*cos(ab));
            M(m,k) = M(m,k)-V(k)*V(m)*(-gkmf*cos(ab)+bkmf*sin(ab));
            HESSIAN_M(m,k)=HESSIAN_M(m,k)-V(k)*V(m)*(-gkmf*sin(ab)-bkmf*cos(ab));
            %%%%%%%%%%%%%%
            L(k,m) = L(k,m)+V(k)*(-gkmf*sin(ab)+bkmf*cos(ab));
            HESSIAN_L(k,m)=HESSIAN_L(k,m)+(-gkmf*sin(ab)+bkmf*cos(ab));
            L(m,k) = L(m,k)-V(m)*(-gkmf*sin(ab)-bkmf*cos(ab));  
            HESSIAN_L(m,k)=HESSIAN_L(m,k)-(-gkmf*sin(ab)-bkmf*cos(ab));
        end

        %Jacobiana
        Jaux = [H N ; M L];
        HESSIAN_aux=[HESSIAN_H HESSIAN_N;HESSIAN_M HESSIAN_L];

        for k=1:nb 
            if tipo(k)==3 
                H(k,k) = 10^10;
                HESSIAN_H(k,k)=10^10;
            end

            if tipo(k)>1; % Barra de geração VTheta ou PV
                L(k,k) = 10^10;
                HESSIAN_L(k,k)=10^10;
            end    
        end    

        J = [H N ; M L];

    %--------------------------------------------------------------------------
    %------------------ Para Criação da Matriz Algébrica ----------------------
    %--------------------------------------------------------------------------

  %  clear H1 N1 M1 L1 Jaux

    H1 = H; N1 = N; M1 = M; L1 = L;

    Jaux = [H1 N1;M1 L1];

    %------------------ Fim da Construção da Matriz Jacobiana -----------------

    %--------------------------------------------------------------------------
    %---------------- Número Infinito para a submatriz H,N,M,L ----------------
    %--------------------------------------------------------------------------

        for i=1:nb
            if (tipo(i) == 2)
                N(:,i) = 0;
                HESSIAN_N(:,i)=0;
                M(i,:) = 0;
                HESSIAN_M(i,:)=0;
                L(i,:) = 0;
                HESSIAN_L(i,:)=0;
                L(:,i) = 0;
                HESSIAN_L(:,i)=0;
                L(i,i) = 1;
                HESSIAN_L(i,i)=1;
            elseif (tipo(i) == 3)
                H(i,:) = 0;
                HESSIAN_H(i,:)=0;
                H(:,i) = 0;
                HESSIAN_H(:,i)=0;
                H(i,i) = 1;     
                HESSIAN_H(i,i)=1;
                N(i,:) = 0;
                HESSIAN_N(i,:)=0;
                N(:,i) = 0;
                HESSIAN_N(:,i)=0;
                M(i,:) = 0;
                HESSIAN_M(i,:)=0;
                M(:,i) = 0;
                HESSIAN_H(:,i)=0;
                L(i,:) = 0;
                HESSIAN_L(i,:)=0;
                L(:,i) = 0;
                HESSIAN_L(:,i)=0;
                L(i,i) = 1;
                HESSIAN_L(i,i)=1;
            end
        end

    %--------------------------------------------------------------------------
    %-------------------- Construção da Matriz Jacobiana ----------------------
    %--------------------------------------------------------------------------

    J=[H N;M L];
    HESSIAN_sup=[HESSIAN_H HESSIAN_N;HESSIAN_M HESSIAN_L];
    %----------------- Fim da Construção da Matriz Jacobiana ------------------
    JAC=J;
    %--------------------------------------------------------------------------
    %------------------------- Resolvendo Jacobiana ---------------------------
    %--------------------------------------------------------------------------

    D = inv(J)*DPQ;

    %--------------------------------------------------------------------------

    %--------------------------------------------------------------------------
    %----------- Obtenção do vetor de correção Delta Teta e Delta V -----------
    %--------------------------------------------------------------------------

    DT = (D(1:nb));  

    DV = (D(nb+1:2*nb)); 

    %--------------------------------------------------------------------------
    %---------------------------- Nova solução --------------------------------
    %--------------------------------------------------------------------------

    T = T + DT; 

    V = V + DV;

    %--------------------------------------------------------------------------
    %------------------- Limpa as submatrizes H M N L -------------------------
    %--------------------------------------------------------------------------

    clear H M N L

    %--------------------------------------------------------------------------
    %----------------- Incrementando Contador de Iterações --------------------
    %--------------------------------------------------------------------------

    if flag == 0    
      it = it + 1; %contador de iterações   
    end

    %--------------------------------------------------------------------------
    %-------------- Rotina para Mostrar o Caminho dos Mismatches --------------
    %--------------------------------------------------------------------------

      iter=[1:it]; %Construção do vetor iteração
      maxdp(it) = max(abs(DP));
      maxdq(it) = max(abs(DQ));

    end % Fim do while

    %--------------------------------------------------------------------------
    %----------------- Gráfico de Mismatches de Pk e QK -----------------------
    %--------------------------------------------------------------------------
%     if plotgraf ==1
%     figure;
%     semilogy(iter,maxdp,'ro -','LineWidth',2);
%     hold on
%     grid on
%     semilogy(iter,maxdq,'k* -','LineWidth',2);
%     %legend = legend('max|Delta P_k|','max |Delta Q_k|');
%     title('Caminho dos mismatches no processo iterativo','FontSize',14);
%     xlabel('Iterações','FontSize',15);
%     ylabel('Tolerância','FontSize',15);
%     end
%     %--------------------------------------------------------------------------
    %--------- Cálculo dos Fluxos e Perdas nas Linhas de Transmissão ----------
    %--------------------------------------------------------------------------

    for i=1:nr
        k=de(i);
        m=para(i);

        %----------------------------------------------------------------------
        akk(i) = (1/(tap(i)^2)); % tap lado de baixa
        akm(i) = 1/tap(i);% tap lado de baixa
        amm(i) = 1.0; 
    end
    % 
  
    for l=1:nr
        k  = de(l);
        m  = para(l);
        ab = T(k) - T(m);

        gkm = r(l)/(r(l)^2+x(l)^2);
        bkm = -x(l)/(r(l)^2+x(l)^2);

        vkm = V(k)*V(m); 
        t   = tap(l);  

        pkm(l) = t^2*V(k)*V(k)*gkm-t*vkm*(gkm*cos(ab)+bkm*sin(ab));
        qkm(l) = -t^2*V(k)*V(k)*(bkm+bshl(l))-t*vkm*(gkm*sin(ab)-bkm*cos(ab));

        pmk(l) = V(m)*V(m)*gkm-t*vkm*(gkm*cos(ab)-bkm*sin(ab));
        qmk(l) = -V(m)*V(m)*(bkm+bshl(l))+t*vkm*(gkm*sin(ab)+bkm*cos(ab));

        pperdas(l) = pkm(l)+pmk(l); 
        qperdas(l) = qkm(l)+qmk(l);  
    end

    %-------- Fim do Cálculo dos Fluxos e Perdas na Linha de Transmissão ------

    %--------------------------------------------------------------------------
    % -------------------- Transformação em MW e MVAR -------------------------
    %--------------------------------------------------------------------------

    %------------------ Fluxos de Potência Ativa e Reativa --------------------

      PkmMW = pkm*MVAbase;
      QkmMVAR = qkm*MVAbase;
      PmkMW = pmk*MVAbase;
      QmkMVAR = qmk*MVAbase;

    %--------------------- Perdas nas Linhas de Trasmissão --------------------

      PeratMW = pperdas*MVAbase;
      PerreatMVAR = qperdas*MVAbase;
      Tpa = sum(PeratMW); %total perdas ativas
      Tpr = sum(PerreatMVAR);%total perdas reativas

    %--------------------- Potências Geradas e Consumidas ---------------------

      PcMW = pc*MVAbase;
      QcMVAR = qc*MVAbase;
      PgMW = (Pcal + pc)*MVAbase;
      QgMVAR = (Qcal + qc)*MVAbase;

%       for i=1:nr  
%           if(i==bi || i==bj || i==bk)
%            PgMW(i) = 0;
%            QgMVAR(i) =0;  
%           end
%       end
      Tpc = sum(PcMW);
      Tqc = sum(QcMVAR);
      Tpg = sum(PgMW);
      Tqg = sum(QgMVAR);

    bus = [(1:nb)',tipobus,PcMW,QcMVAR,gshk,bshk,narea,V,T*radgraus];

    save bus

    branch = [de, para, r, x, bshl*2];

    save branch

    ybarra = Y;
    ybarra_semshuntBL=YY;

    save ybarra
    save ybarra_semshuntBL
    
    ngaux=1;
    ibarra_auxiliar=zeros();
    PgMW_auxiliar=zeros();
    QgMVAR_auxiliar=zeros();
    for ibarra=1:nb
        if tipobus(ibarra)==2|| tipobus(ibarra)==3 
        ibarra_auxiliar(ngaux)=ibarra;
        PgMW_auxiliar(ngaux)=PgMW(ibarra);
        QgMVAR_auxiliar(ngaux)=QgMVAR(ibarra);
            ngaux=ngaux+1;
        end
    end
    ng=ngaux;
    %ng = length(find(pg));

    gen = [ibarra_auxiliar' PgMW_auxiliar' QgMVAR_auxiliar'];

    save gen
    save V

    %------------------- Fim da Transformação em MW e MVAR --------------------  

    %--------------------------------------------------------------------------
    %---------------------- Gerando Relatório Final ---------------------------
    %--------------------------------------------------------------------------
% 
%     c = fopen('Relatorio.out','w+');
%     format short 
%     fprintf(c,'\n O Fluxo de Carga Pelo Método de Newton-Raphson convirgiu em %g iterações\n\n',it);
%     fprintf(c,'====================================================================================\n');
%     fprintf(c,'-------------------------------- SUMÁRIO DO SISTEMA --------------------------------\n');
%     fprintf(c,'====================================================================================\n');
%     fprintf(c,'           Quantidades?                 |              Outras Informações           \n');
%     fprintf(c,'------------------------------------------------------------------------------------\n');
%     fprintf(c,'     Barras:             %1.0f',nb);fprintf(c,'                   Tensão Máxima:             %1.3f \n',(max(V)));            
%     fprintf(c,'     Ramos:              %1.0f',nr);fprintf(c,'                   Tensão Mínima:             %1.3f \n',(min(V)));  
%     fprintf(c,'     Geradores:           %1.0f',(nslack+npv));fprintf(c,'                   Ângulo Máximo:             %1.3f \n',max((T*radgraus))); 
%     fprintf(c,'     Cargas:             %1.0f',ncargas);fprintf(c,'                   Ângulo Mínimo:            %1.3f \n',min((T*radgraus)));
%     fprintf(c,'     Trafos:              %1.0f',ntransf);fprintf(c,'                   P Perdas (I^2*R):          %1.3f \n',Tpa);
%     fprintf(c,'     Tolerância:         %1.g',tol);fprintf(c,'               Q Perdas (I^2*X):          %1.3f \n',Tpr); 
%     fprintf(c,'====================================================================================\n');
% 
%     fprintf(c,'====================================================================================\n\n\n');
% 
%     %--------------------------------------------------------------------------
%     %-------------------- Tensões e Ângulos das Barras ------------------------
%     %--------------------------------------------------------------------------
%     fprintf(c,'====================================================================================\n');
%     fprintf(c,'------------------------------ RELATÓRIO DAS TENSÕES -------------------------------\n');
%     fprintf(c,'====================================================================================\n');
%     fprintf(c,'| Barra | Tipo |    Mag    |  Fase(°)  | PG(MW)  | QG(MVAr)|  PC(MW) | QC(MVAr)|\n');
%     fprintf(c,'====================================================================================\n');
%     format short  
%     for k = 1:nb
%         fprintf(c,'%7d | %4d | %9.4f | %9.4f | %7.2f | %7.2f | %7.2f | %7.2f | \n',...
%             numb(k),tipo(k),V(k),(T(k)*radgraus),PgMW(k),QgMVAR(k),PcMW(k),QcMVAR(k));
%     end
%     fprintf(c,'------------------------------------------------------------------------------------\n');
%     fprintf(c,'|                               Total: | %7.2f | %7.2f | %7.2f | %7.2f |\n',...
%             Tpg,Tqg,Tpc,Tqc);
%     fprintf(c,'====================================================================================\n');
%     
%     fprintf(c,'====================================================================================\n\n\n');
% 
% 
%     %--------------------------------------------------------------------------
%     % ------------------- Impressão dos Fluxos nas Linhas ---------------------
%     %--------------------------------------------------------------------------
% 
%     fprintf (c,'====================================================================================\n');
%     fprintf (c,'-------------------------- RELATÓRIO DE FLUXOS NAS LINHAS --------------------------\n');
%     fprintf (c,'------------------------------------------------------------------------------------\n');
%     fprintf (c,'|  DE-PARA |   Pkm(MW) |  Qkm(MW)  | Pmk(MVAR) | Qmk(MVAR) |   PA(MW) |  PR(MVAR)  |\n');
%     fprintf (c,'------------------------------------------------------------------------------------\n');
%     format short
%     for i=1:nr    
%     %if(para(i)~=bj && para(i)~=bk)    
%     fprintf(c,'| %2g -%3g  |%9.2f  |%9.2f  |%9.2f  |%9.2f  |%8.3f   | %8.2f  |',...
%     de(i),para(i),PkmMW(i),QkmMVAR(i),PmkMW(i),QmkMVAR(i),PeratMW(i),PerreatMVAR(i));
%     fprintf(c,'\n');  
%     %end
%     end
%     fprintf(c,'------------------------------------------------------------------------------------\n');
%     fprintf(c,'|                                                    Total:  %7.3f      %7.2f  |\n',Tpa,Tpr);
%     fprintf(c,'====================================================================================\n');
%    
%     fprintf(c,'====================================================================================\n');
% 
%     open Relatorio.out
 end