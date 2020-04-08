function AstroSeis_liquidcore(varargin)
% program to compute seismic wavefield in a solid body with a liquid core  
% Yuan Tian 03/17/2020 @University of Houston
if numel(varargin)<1
    fprintf('AstroSeis_liquidcore need input parameter file !\n');
    return;
end

set(0,'defaultaxesfontsize',14);
set(0,'defaulttextfontsize',14);
set(0,'defaultLineLineWidth',1);
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
        load(mesh_file_name);
        tmpline = fgetl(tmpfid);
        tmpline = fgetl(tmpfid);
        tmpline = fgetl(tmpfid);
    else
        tmpline = fgetl(tmpfid);
        tmp = textscan(tmpline,'%f %f %f','CommentStyle','#');
        rpl(1)=tmp{1};
        nmesh(1)=tmp{2};
        nfold(1)=tmp{3};
        tmpline = fgetl(tmpfid);
        tmp = textscan(tmpline,'%f %f %f','CommentStyle','#');
        rpl(2)=tmp{1};
        nmesh(2)=tmp{2};
        nfold(2)=tmp{3};
        fprintf('generating mesh \n');
        fprintf('body radius=%f, nmesh=%d \n',rpl(1),nmesh(1));
        [layer]=gen_layer(rpl,nmesh,nfold);
        tmpline = fgetl(tmpfid);
        tmp = textscan(tmpline,'%s','CommentStyle','#');
        out_mesh_name=tmp{1}{1};
        save(out_mesh_name,'layer');
    end
    tmpline = fgetl(tmpfid);
    tmp = textscan(tmpline,'%s','CommentStyle','#');
    output_file_name=tmp{1}{1};
    tmpline = fgetl(tmpfid);
    tmp = textscan(tmpline,'%f %f %f %f','CommentStyle','#');
    vp1=tmp{1};
    vs1=tmp{2};
    rho1=tmp{3};
    Q=tmp{4};
    tmpline = fgetl(tmpfid);
    tmp = textscan(tmpline,'%f %f %f %f','CommentStyle','#');
    vp2=tmp{1};
    vs2=tmp{2};
    rho2=tmp{3};
    tmpline = fgetl(tmpfid);
    tmp = textscan(tmpline,'%f %f %f','CommentStyle','#');
    nt=tmp{1};
    dt=tmp{2};
    f0=tmp{3};
    T=nt*dt; % total length of time
    fmax=3*f0; % maxium frequency to compute
    disp(['Pwavelength ',num2str(vp1/f0),'(m)']);
    disp(['Swavelength ',num2str(vs1/f0),'(m)']);
    disp(['min_ele_size',num2str(min([layer(2).face.a])),'(m)']);
    disp(['max_ele_size',num2str(max([layer(2).face.a])),'(m)']);
    tmpline = fgetl(tmpfid);
    tmp = textscan(tmpline,'%s %f','CommentStyle','#');
    source_type=tmp{1}{1};
    source_scale=tmp{2};
    if strcmp(source_type,'single')
        tmpline = fgetl(tmpfid);
        tmp = textscan(tmpline,'%f %f %f','CommentStyle','#');
        fsrc=10^source_scale*[tmp{1},tmp{2},tmp{3}];
        tmpline = fgetl(tmpfid);
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
 face1=layer(2).face;
 face2=layer(1).face;
for i=1:length(face1)
    Ra(i)=mean(norm(face1(i).ic));
end
for i=1:length(face2)
    Rca(i)=mean(norm(face2(i).ic));
end
R=mean(Ra);
rs = R-h; 
Rc=mean(Rca);
xs=rs*sin(thetas)*cos(phis);
ys=rs*sin(thetas)*sin(phis);
zs=rs*cos(thetas);
wi=4/T; %imaginary part angular frequency
iu=complex(0,1);
df = 1/(nt*dt);
mu1 =rho1*vs1*vs1;
lamda1 = rho1*vp1*vp1 - 2*mu1;
mu2 =rho1*vs2*vs2;
lamda2 = rho2*vp2*vp2 - 2*mu2;
nint = 10;
nxi=300;
%% transform matrix
Smat=Smat_func(face2);
%% computing

for iw = 2: nt
iw
    tic
    freq = (iw-1)*df
    if freq > fmax
        continue;
    end
    %wi=1i*pi/T;
    w = 2*pi*freq+wi*iu;
    if strcmp(source_type,'single')
        if rs<Rc
            disp('single force source is in liquid part, change source location')
            break;
        else
            u02=u0e_func(face2,w,rho1,mu1,lamda1,xs,ys,zs,Q,fsrc);
            u01=u0e_func(face1,w,rho1,mu1,lamda1,xs,ys,zs,Q,fsrc);
            P0=zeros(length(face2),1);
        end
    elseif strcmp(source_type,'moment')
        if rs<Rc
            u02=zeros(length(face2)*3,1);
            u01=zeros(length(face1)*3,1);
            P0=u0p_func_exp(face2,w,vp2,xs,ys,zs,Q);
        else
            u02=u0eM_func(face2,w,rho1,mu1,lamda1,xs,ys,zs,Q,M);
            u01=u0eM_func(face1,w,rho1,mu1,lamda1,xs,ys,zs,Q,M);
            P0=zeros(length(face2),1);
        end
    end
    [A,b]=liq_core(face1,face2,w,lamda1,mu1,rho1,lamda2,mu2,rho2,Q,nint,...
        nxi,P0,u01,u02,Smat,2*pi*f0);
    x=A\b;
    up(:,iw)=x(1:length(face2));
    u2(:,iw)=x(length(face2)+1:4*length(face2));
    u1(:,iw)=x(4*length(face2)+1:end);
    %uuu(:,iw)=u02;
    %uu(:,iw)=P0;
    if iw >= 2 
        up(:,nt+2-iw) = conj(up(:,iw)); 
        u2(:,nt+2-iw) = conj(u2(:,iw)); 
        u1(:,nt+2-iw) = conj(u1(:,iw)); 
    end
    save(output_file_name,'up','u2','u1','nt','T');
    
     toc
end


end