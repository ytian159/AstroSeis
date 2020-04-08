function [layer]=gen_layer(rpl,nmesh,nfold)
nlayer=length(rpl);
transpar=0.7;
%%
figure
hold on;
for i=1:nlayer
    [face,numface,ds,xs0,ys0,zs0,thetas,phis,h,V,Tri]=gen_mesh_topo_layers(rpl(i),nmesh(i),nfold(i),transpar);
    layer(i).face=face;
    layer(i).numface=numface;
    layer(i).ds=ds;
    layer(i).xs0=xs0;
    layer(i).ys0=ys0;
    layer(i).zs0=zs0;
    layer(i).thetas=thetas;
    layer(i).phis=phis;
    layer(i).h=h;
    layer(i).V=V;
    layer(i).Tri=Tri;
    
end
xlabel('X (m)')
ylabel('Y (m)')
zlabel('Z (m)')
pause(0.1);
hold off
end