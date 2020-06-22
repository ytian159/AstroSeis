function [ x, y, w ] = rule20 ( )

%*****************************************************************************80
%
%% RULE20 returns the rule of degree 20.
%
%  Discussion:
%
%    Order 20 (79 pts)
%    1/6 data for 20-th order quadrature with 18 nodes.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    26 June 2014
%
%  Author:
%
%    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
%    This MATLAB version by John Burkardt.
%
%  Parameters:
%
%    Output, real X(*), Y(*), the coordinates of the nodes.
%
%    Output, real W(*), the weights.
%
  x = [ ...
      0.00000000000000000000000000000000, ...
      0.00000000000000000000000000000000, ...
      0.00000000000000000000000000000000, ...
     -0.86696389117550812661190351318543, ...
     -0.46255868049954115157854976204512, ...
      0.00000000000000000000000000000000, ...
     -0.67416180418116902260753583861473, ...
     -0.22447168033665606945391317277106, ...
      0.00000000000000000000000000000000, ...
      0.00000000000000000000000000000000, ...
     -0.55640337063475932672618703937903, ...
      0.00000000000000000000000000000000, ...
     -0.76173172264834809386157300085862, ...
     -0.15012093407474928012860293037678, ...
     -0.27874288623783821051052255615910, ...
     -0.42809996429665844299450237790398, ...
      0.00000000000000000000000000000000, ...
      0.00000000000000000000000000000000 ];
  y = [ ... ...
      0.50935573580030290124442500215952, ...
      0.10254518566344431484070667698816E+01, ...
     -0.49506265376045875810405496904400, ...
     -0.56894127058564452950719739998646, ...
     -0.39335935346809758365690922062188, ...
     -0.38873359764810130895317005413822, ...
     -0.56423729270253943031731484496373, ...
     -0.33519558519158101515864503661923, ...
      0.27281208605145075547218518170954, ...
      0.00000000000000000000000000000000, ...
     -0.49670535155060417023204843137510, ...
     -0.20816484443009696237324462486882, ...
     -0.51090241799312035370196362374051, ...
     -0.56032152802942467258347810157577, ...
     -0.48210916153383866244702005141785, ...
     -0.55875287099133282930781054909746, ...
      0.11166780705147990599713407600049E+01, ...
      0.77578464434062421187747427697944 ];
  w = [ ... ...
      0.12072956229196140184720818378917E-01, ...
      0.28443984028101927395806965629882E-02, ...
      0.93465277263442986124162626755397E-02, ...
      0.29739840427656477081802036308844E-02, ...
      0.20327046933776882823607497409423E-01, ...
      0.12440057912161099228836551604472E-01, ...
      0.57983520915299366578479709959397E-02, ...
      0.30774405447413411002155380888934E-01, ...
      0.18534535260005961010495899016617E-01, ...
      0.61022450704916040183649592581025E-02, ...
      0.15757087201875994349887096272192E-01, ...
      0.18146095122192897098014042659411E-01, ...
      0.10912126413380354963231665469691E-01, ...
      0.97275807652705556020653743902143E-02, ...
      0.22813420666829580712213934309259E-01, ...
      0.94183526939491428417826897656981E-02, ...
      0.10513336056091899919123590224643E-02, ...
      0.10305163239812520591223081322085E-01 ];

  return
end
