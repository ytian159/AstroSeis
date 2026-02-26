function Mesh_generater_lc(varargin)
% program to compute mesh for a model in a solid body with a liquid core  
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
    
    % input mesh file
    tmpline = fgetl(tmpfid);
    mesh_file_name =  textscan(tmpline,'%s','CommentStyle','#');
    mesh_file_name = char(mesh_file_name{1});
    if exist(mesh_file_name,'file')
        disp('the mesh file already exist \n');
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
end
end