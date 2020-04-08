function wiggler(x, y, data, sc_fac, pos_neg,nfold)
%FUNCTION WIGGLE(X,Y,DATA, SC_FAC, POS_NEG) creates a wiggle plot (e.g. for
%seismic or radar data), where either the 'positive' or 'negative' part of
%the wavelets is filled.
%Application: - WIGGLE(DATA) plots data by assuming X=1:NX and Y=1:NY,
%               where [NY,NX]=SIZE(DATA), SC_FAC=1 and POS_NEG='P'
%             - WIGGLE(DATA,SC_FAC) or WIGGLE(DATA,POS_NEG) either assumes
%               POS_NEG = 'P' or SC_FAC = '1'. The remaining unknows are
%               computed as above.
%             - WIGGLE(DATA,SC_FAC,POS_NEG) plots data by assuming
%               X=1:NX and Y=1:NY (as above).
%             - WIGGLE(X,Y,DATA) plots wiggles filled above zero
%               (pos_neg='p')
%             - WIGGLE(X,Y,DATA,POS_NEG) or WIGGLE(X,Y,DATA,SC_FAC), the
%               remaining unknown is assumed accoring above examples.
%
%IMPORTANT: the function assumes that the amount of samples along the
%'time' axis (y-axis) is larger then the amount of traces to plot! It is
%thus possible to control the orientation of the traces by transposing the
%data-matrix.
%E.g. x=1:10; y=1:1000; data=randn(1000,10); -> traces are plotted along
%     x-axis (x-axis==time-axis).
%
%REMARK: it lies in the nature of this function that it plots a white line
%to cover interpolation effects of the fill-function. It is therefore not
%advisable to overplot an image (e.g. created with imagesc) with
%wiggled-traces, as a white line would show up! To avoid this drawback, the
%only solution (I can think of, at the moment) would be to use the
%patch-function instead of the fill-function, rendering the wiggle-creation
%to a very time consuming process...
%
%(c) Jacques Ernst 2006                                         ETH Zurich
%=========================================================================

if (nargin==1)
    data = x; clear x;
    [ny,nx] = size(data); x = 1:nx; y=1:ny; sc_fac=1; pos_neg='p';
end
if (nargin==2)
    data = x; clear x;
    if (isstr(y)) pos_neg = y; clear y; sc_fac = 1;
    else sc_fac = y; clear y; pos_neg = 'p'; end
    [ny,nx] = size(data); x = 1:nx; y=1:ny;
end
if (nargin==3 && isstr(data))
    pos_neg = data; clear data;
    sc_fac = y; clear y;
    data = x; clear x;
    [ny,nx] = size(data); x = 1:nx; y=1:ny;
elseif (nargin==3)
    [ny,nx] = size(data);
    sc_fac = 1; pos_neg = 'p';
end
if (nargin==4 && isstr(sc_fac))
    pos_neg = sc_fac; clear sc_fac; sc_fac = 1;
elseif (nargin==4)
    [ny,nx] = size(data); pos_neg = 'p';
end
if (nargin==5) [ny,nx] = size(data);
end
if (nargin==6) [ny,nx] = size(data);
end
if (min([nx ny])>1) M_data = max(data(:)); data = data./M_data; end

idx_pos = find(data>=0);
idx_neg = find(data<=0);

switch lower(pos_neg)
    case { 'p' }
        data_fill = data; data_fill(idx_neg) = 0.0;
    case { 'n' }
        data_fill = data; data_fill(idx_pos) = 0.0;
end
data=data*nfold;
%Choose plot-direction depending on the nx,ny-sizes (automatically)
if (nx<ny)
    %first & last points in data (for each trace) must be 0.0 in order to
    %create a clean filling!
    data_fill(:,1) = 0.0; data_fill(length(data),:) = 0.0;
    for (ii=1:nx)
        %handle_fill = fill(data_fill(:,ii).*sc_fac+x(ii)-1,y,'k');
        %set(handle_fill,'LineStyle','none');
        hold on;
        %plot(ones(1,length(y)).*x(ii)-1,y,'w',data(:,ii).*sc_fac+x(ii),y,'r');
        plot(data(:,ii).*sc_fac+x(ii),y,'r');
    end
else
    %last point in data (for each trace) must be 0.0 in order to create a clean
    %filling!
    data_fill(:,1) = 0.0; data_fill(:,length(data)) = 0.0;
    for (ii=1:ny)
        %handle_fill = fill(x,data_fill(ii,:).*sc_fac+y(ii)-1,'k');
        %set(handle_fill,'LineStyle','none');
        hold on
        %plot(x,ones(1,length(x)).*y(ii)-1,'w',x,data(ii,:).*sc_fac+y(ii),'r');
        plot(x,data(ii,:).*sc_fac+y(ii),'r');
    end
end
hold off;
axis tight;

