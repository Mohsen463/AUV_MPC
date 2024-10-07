
%% // Default parameters
param.meshsize  = 256 ;     %// main grid size
param.patchsize = 100 ;     
param.windSpeed = 50  ;    %// what unit ? [m/s] ??
param.winddir   = 90   ;    %// Azimuth
param.rng       = 05 ;            %// setting seed for random numbers
param.A         = 01 ;    %// Scaling factor
param.g         = 9.81 ;    %// gravitational constant

%% // get the grid parameters which remain constants (not time dependent)
[H0, W, Grid_Sign] =  initialize_wave( param ) ;
H1 = (H0*1e4);

%% // calculate wave at t0
t0 = 0 ;
Z0 = calc_wave(  H1 , W , t0 , Grid_Sign ) ;

 % figure(1), mesh(Z0)
 % figure(2), mesh(Z1*100)

for ii = 1: 256 
    for jj = 4:4: 256*4
        Z(ii, jj-3:jj) = Z0(ii,jj/4);
    end
end
Z1=Z;
for ii = 6:6:256*6
    for jj = 1:256*4
        Z(ii-5:ii, jj) = Z1(ii/6,jj);
    end
end
Z=Z-depth+100;

Z(1:1536,1025:1024+300) = Z(1:1536,1:300);
Z=Z';
% figure(3), mesh(Z)
% zlim([-700 0])

% Z2 (1:1000,1:1000) = 0;
% Z2(245:244+512,245:244+512)=Z;
% figure(4), mesh(Z2)

%% // populate the display panel
% h.fig  = figure('Color','w') ;
% h.ax   = handle(axes) ;                 %// create an empty axes that fills the figure
% h.surf = handle( surf( NaN(2) ) ) ;     %// create an empty "surface" object
% 
% %% // Display the initial wave surface
% set( h.surf , 'XData',X , 'YData',Y , 'ZData',Z )
% set( h.ax   , 'XLim',param.xLim , 'YLim',param.yLim , 'ZLim',param.zLim )
% 
% %% // Change some rendering options
% % axis off                                %// make the axis grid and border invisible
% shading interp                          %// improve shading (remove "faceted" effect)
% blue = linspace(0.4, 1.0, 25).' ;% cmap = [blue*0, blue*0, blue]; %'// create blue colormap
% % colormap(cmap)
% %// configure lighting
% h.light_handle = lightangle(-45,30) ;   %// add a light source
% set(h.surf,'FaceLighting','phong','AmbientStrength',.3,'DiffuseStrength',.8,'SpecularStrength',.9,'SpecularExponent',25,'BackFaceLighting','unlit')

%% // Animate
% view(75,55) %// no need to reset the view inside the loop ;)
% 
% timeStep = 1./25 ;
% nSteps = 2000 ;
% for time = (1:nSteps)*timeStep    
%     %// update wave surface
%     Z = calc_wave( H0,W,time,Grid_Sign ) ;
%     h.surf.ZData = Z ;
%     pause(0.001);
% end


%% // This block of code is only if you want to generate a GIF file
%// be carefull on how many frames you put there, the size of the GIF can
%// quickly grow out of proportion ;)
% 
% nFrame = 55 ;
% gifFileName = 'MyDancingWaves.gif' ;
% 
% view(-70,40)
% clear im
% f = getframe;
% [im,map] = rgb2ind(f.cdata,256,'nodither');
% im(1,1,1,20) = 0;
% iframe = 0 ;
% for time = (1:nFrame)*.5
%     %// update wave surface
%     Z = calc_wave( H0,W,time,Grid_Sign ) ;
%     h.surf.ZData = Z ;
%     pause(0.001);
% 
%     f = getframe;
%     iframe= iframe+1 ;
%     im(:,:,1,iframe) = rgb2ind(f.cdata,map,'nodither');
% end
% imwrite(im,map,gifFileName,'DelayTime',0,'LoopCount',inf)
% disp([num2str(nFrame) ' frames written in file: ' gifFileName])