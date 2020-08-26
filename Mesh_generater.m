function Mesh_generater(varargin)
% program to generate mesh
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
        disp('the mesh file already exist \n');
        tmpline = fgetl(tmpfid);
        tmpline = fgetl(tmpfid);
    else
        tmpline = fgetl(tmpfid);
        tmp = textscan(tmpline,'%f %f %f %s','CommentStyle','#');
        R=tmp{1};
        nmesh=tmp{2};
        nfold=tmp{3};
        aster_name=tmp{4};
        fprintf('generating mesh \n');
        fprintf('R=%f, nmesh=%d \n',R,nmesh);
        if strcmp(aster_name,'rand')
            [face,numface,ds,xs0,ys0,zs0,thetas,phis,height,V,Tri]=gen_mesh_topo_rand(R,nmesh,nfold);
        else
            [face,numface,ds,xs0,ys0,zs0,thetas,phis,height,V,Tri]=gen_mesh_ph_topo(R,nmesh,nfold);
        end
        tmpline = fgetl(tmpfid);
        tmp = textscan(tmpline,'%s','CommentStyle','#');
        out_mesh_name=tmp{1}{1};
        save(out_mesh_name,'face','numface','ds','xs0','ys0','zs0','thetas','phis','height','V','Tri');
    end
    
end
end