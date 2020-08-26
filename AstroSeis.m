function AstroSeis(varargin)
% program to compute seismic wavefield in homogenous elastic medium
% Yuan Tian 03/15/2020 @University of Houston
if numel(varargin)<1
    fprintf('AstroSeis need input parameter file !\n');
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
    %% plot mesh
    for i=1:length(face)
        Ra(i)=mean(face(i).ic);
    end
    R=mean(Ra);
    nV=length(V);
    for i=1:nV
        hh=norm(V(i,:));
    end
    figure; hold on;
    h= patch('faces',Tri,'vertices',V, 'FaceVertexCData', hh(:), 'FaceColor','interp');
    colormap (jet);
    alpha(h,.6);
    set(h,'EdgeColor',[152 57 153]/255,'linewidth',.01);
    % set(h,'EdgeColor','b','FaceColor',[1 1 1 ]*.5)
    axis equal vis3d
    view(3)
    colorbar;
    pause(.1);
    %%
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

rs = R-h; 
xs=rs*sin(thetas)*cos(phis);
ys=rs*sin(thetas)*sin(phis);
zs=rs*cos(thetas);
wi=4/T;
iu=complex(0,1);
df = 1/(nt*dt);
mu =rho*vs*vs;
lamda = rho*vp*vp - 2*mu;
%% computing
uu=zeros(length(face)*3,nt);
for iw = 2: nt
iw
    tic
    freq = (iw-1)*df
    if freq > fmax
        break;
    end
    w = 2*pi*freq+wi*iu;
    if strcmp(source_type,'single')
        [u01]=u0e_func(face,w,rho,mu,lamda,xs,ys,zs,Q,fsrc);
    elseif strcmp(source_type,'moment')
        [u01]=u0eM_func(face,w,rho,mu,lamda,xs,ys,zs,Q,M);
    end
    TRAC=cal_traction_tri_vec(face,length(face),w,lamda,mu,rho,Q);
    uu(:,iw)=TRAC\u01;
    if iw >= 2 
        %uu(:,nt+2-iw) = conj(uu(:,iw)); 
        uu(:,nt+2-iw) = conj(uu(:,iw)); 
    end
     save(output_file_name,'uu','nt','T');
     toc
end


end