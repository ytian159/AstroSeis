# AstroSeis
his MATLAB codes compute 3-D seismic wavefield using BEM.
For homogenous model, run:

AstroSeis BEM_para

in MATLAB. For solid body with a liquid core model, run:

AstroSeis_liquidcore  BEM_para_lc

in MATLAB.Visualization is also embedded in this code. 
The input file can be the same with  the computation part. 
If the computation is finished, we can simply run 

AstroSeis_plot    Phobos_topo_model 

to do the plotting of the homogenous model. 
For the liquid core model is similar: 

AstroSeis_liquidcore_plot   inp_sft_lc

There are three examples that is already finished in "exmaples" folder. 
Usage is similar, you only need to change the input files. 
If you have any questions, please feel free to contact to:
ytian159@gmail.com
Yuan Tian
