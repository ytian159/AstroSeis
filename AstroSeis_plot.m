function AstroSeis_plot(varargin)
% program to plot seismic wavefield in homogenous elastic medium
% Yuan Tian 03/25/2020 @University of Houston
if numel(varargin)<1
    fprintf('AstroSeis need input parameter file !\n');
    return;
end
set(0,'defaultaxesfontsize',25); % control the default fontsize of the plots
set(0,'defaulttextfontsize',25);
set(0,'defaultlinelinewidth',1.5);
%% set up path
pathAS=fileparts(mfilename('fullpath'));
if contains(pathAS,'examples')
    pathAS=extractBefore(path,'examples');
end
addpath(genpath(pathAS));
%addpath('lib_BEM');
parafilefn = varargin{1};

%% read in the parameter file 
if ~exist(parafilefn,'file')
    fprintf('%s  is not exist ! \n',parafilefn);
else
    tmpfid = fopen(parafilefn,'r');
    
    % input segyfile 
    tmpline = fgetl(tmpfid);
    mesh_file_name =  textscan(tmpline,'%s','CommentStyle','#');
    mesh_file_name = char(mesh_file_name{1});
    if exist(mesh_file_name,'file')
        load(mesh_file_name,'face','Tri','V');
        tmpline = fgetl(tmpfid);
        tmpline = fgetl(tmpfid);
    else
        tmpline = fgetl(tmpfid);
        tmp = textscan(tmpline,'%f %f %f','CommentStyle','#');
        R=tmp{1};
        nmesh=tmp{2};
        nfold=tmp{3};
        fprintf('generating mesh \n');
        fprintf('R=%f, nmesh=%d \n',R,nmesh);
        [face,numface,ds,xs0,ys0,zs0,thetas,phis,height,V,Tri]=gen_mesh_ph_topo(R,nmesh,nfold);
        tmpline = fgetl(tmpfid);
        tmp = textscan(tmpline,'%s','CommentStyle','#');
        out_mesh_name=tmp{1}{1};
        save(out_mesh_name,'face','numface','ds','xs0','ys0','zs0','thetas','phis','height','V','Tri');
    end
    tmpline = fgetl(tmpfid);
    tmp = textscan(tmpline,'%s','CommentStyle','#');
    output_file_name=tmp{1}{1};
    tmpline = fgetl(tmpfid);
    tmp = textscan(tmpline,'%f %f %f %f','CommentStyle','#');
    vp=tmp{1};
    vs=tmp{2};
    rho=tmp{3};
    Q=tmp{4};
    tmpline = fgetl(tmpfid);
    tmp = textscan(tmpline,'%f %f %f','CommentStyle','#');
    nt=tmp{1};
    dt=tmp{2};
    f0=tmp{3};
    T=nt*dt; % total length of time
    fmax=3*f0; % maxium frequency to compute
    disp(['Pwavelength ',num2str(vp/f0),'(m)']);
    disp(['Swavelength ',num2str(vs/f0),'(m)']);
    disp(['min_ele_size',num2str(min([face.a])),'(m)']);
    disp(['max_ele_size',num2str(max([face.a])),'(m)']);
    tmpline = fgetl(tmpfid);
    tmp = textscan(tmpline,'%s %f','CommentStyle','#');
    source_type=tmp{1}{1};
    source_scale=tmp{2};
    if strcmp(source_type,'single')
        tmpline = fgetl(tmpfid);
        tmp = textscan(tmpline,'%f %f %f','CommentStyle','#');
        fsrc=10^source_scale*[tmp{1},tmp{2},tmp{3}];
    elseif strcmp(source_type,'moment')
        tmpline = fgetl(tmpfid);
        tmpline = fgetl(tmpfid);
        tmp = textscan(tmpline,'%f %f %f %f %f %f','CommentStyle','#');
        M=10^source_scale*[tmp{1},tmp{2},tmp{3};
            tmp{2},tmp{4},tmp{5};
            tmp{3},tmp{5},tmp{6}];
    else
        disp('source type not recognized\n');
        return;
    end
    tmpline = fgetl(tmpfid);
    tmp = textscan(tmpline,'%f %f %f','CommentStyle','#');
    h=tmp{1};
    thetas = 90-tmp{2}*pi/180;
    phis = tmp{3}*pi/180;
    
    fclose(tmpfid); 
end

%% other parameters
for i=1:length(face)
    Ra(i)=mean(face(i).ic);
end
R=mean(Ra);
rs = R-h; 
xs=rs*sin(thetas)*cos(phis);
ys=rs*sin(thetas)*sin(phis);
zs=rs*cos(thetas);
wi=4/T;
iu=complex(0,1);
df = 1/(nt*dt);
mu =rho*vs*vs;
lamda = rho*vp*vp - 2*mu;
load(output_file_name,'uu','nt','T');
ntr=72;
time = [0:nt-1]*dt;
df = 1/(nt*dt);
dw = 2*pi*df;
ts = 1.5/f0; 
[w,t] = ricker(f0,ts,dt,nt); 
wf = ifft(w)*dt*nt;
 for i=1:length(face)
    xs0(i)=face(i).ic(1);
    ys0(i)=face(i).ic(2);
    zs0(i)=face(i).ic(3);
    radius(i)=norm(face(i).ic);
 end
ngd1=length(face);
%% plotting
u1=uu.';
%%
for j=1:length(face)*3
    %uu2(:,j)=uu(:,j).*wf(:);
    uu2(:,j)=u1(:,j).*wf(:);
end
%uut2_nl=real(fft(uu2_nl))*dw/(2*pi);
uut2=real(fft(uu2))*dw/(2*pi);


phiar=(1:ntr)/ntr*2*pi;
nts=round(ts/dt);
tplot=nt;
ntpp=length((1+nts:tplot));
%nts=0
ub1=zeros(ntpp,ntr);
ub2=zeros(ntpp,ntr);
ub3=zeros(ntpp,ntr);
loc_grid(:,1)=xs0;
loc_grid(:,2)=ys0;
loc_grid(:,3)=zs0;
% figure
% hold on;
ta=t(1+nts:tplot)-ts;
for i=1:ntr
    % receiver location
    rr = R;
    thetar = 0.5*pi;
    phir =phiar(i);
    xr=rr*sin(thetar)*cos(phir);
    yr=rr*sin(thetar)*sin(phir);
    zr=rr*cos(thetar);
    vecr=[xr;yr;zr];
    normvecr=loc_grid*vecr/rr^2;

    [vdum,ind_rec]=max(normvecr);
    xr=xs0(ind_rec); yr=ys0(ind_rec); zr=zs0(ind_rec);
    %plot3(xr,yr,zr,'*');
    ub1(:,i)=uut2((1+nts:tplot),ind_rec);
    ub2(:,i)=uut2((1+nts:tplot),ind_rec+ngd1);
    ub3(:,i)=uut2((1+nts:tplot),ind_rec+2*ngd1);
    vecrn=vecr/norm(vecr);
    uz(:,i)=ub1(:,i)*vecrn(1)+ub2(:,i)*vecrn(2)+ub3(:,i)*vecrn(3);
    %uz(:,i)=uz(:,i).*exp(ta.'/T);
end

%%

figure
hold on
nfold=9;
wiggle(ta,(1:ntr)*5,uz(:,1:ntr).',1,'p',nfold);
xlim([0 40])
xlabel('time (s)');
ylabel('longtitude (°)')
title('Wiggle plot')

%h1=plot(linspace(60,70,ndata),5*ones(ndata,1),'r-');
%%
figure
imagesc(ta,(1:ntr)*5,uz(:,1:ntr).')
xlim([0 40]);

ylabel('longtitude (°)')
colormap jet
umax=max(abs(uz(:)));
colorbar
caxis([-umax umax]*0.2)
title('Event plot');


end