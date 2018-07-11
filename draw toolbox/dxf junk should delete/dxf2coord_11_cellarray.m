%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %                                                     
% dxf2coord 1.1 cell_array                                                %
% author: lukas wischounig, innsbruck, austria (dept. of geology,         %
% university innsbruck), email: csad0018@uibk.ac.at                       %
% date: may 2005                                                          %
% filename: dxf2coord_11_cellarray.m                                      %
% platform: matlab r14                                                    %
% input: dxf versions r2000 - r2004                                       %
% output: cell arrays for 3dpoly's,polylines and lines + matrix for points%
% output form: id   x-coords   y-coords   z-coords                        %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       see !readme.m for details                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
[filename, pathname] = uigetfile({'*.dxf'},'Open File','Multiselect','off'); % choose file to open
addpath(pathname); % add path to the matlab search path
fid=fopen(filename); % open file
C=textscan(fid,'%s'); % read dxf file as cell array of strings C
fclose(fid); % close file to accelerate further computation
C=C{1,1}; % reshape array

%%%%%%%%%%%%%%%%%%%%%%%%%% case point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get some markers and help variables
indpoint=strcmp('AcDbPoint', C); % get line no. of points
pointnum=sum(indpoint); % get total number of lines
indpoint=find(indpoint == 1); % get line no. of lines
points=cell(pointnum,1); % preallocate variable to increase speed
for i=1:pointnum
    point=zeros(1,4);
    point(1,1)=i; % id of line 
    point(1,2)=str2double(C(indpoint(i)+2)); % x start
    point(1,3)=str2double(C(indpoint(i)+4)); % y start
    point(1,4)=str2double(C(indpoint(i)+6)); % z start  
    points{i}=point;
end
clear indpoint point pointnum  % delete garbage from workspace
%%%%%%%%%%%%%%%%%%%%%%%%%% end case point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% case line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get some markers and help variables
indline=strcmp('AcDbLine', C); % get line no. of LW polylines
linenum=sum(indline); % get total number of lines
indline=find(indline == 1); % get line no. of lines
lines=cell(linenum,1); % preallocate variable to increase speed
for i=1:linenum
    ten=strcmp('10',C(indline(i):indline(i)+2));
    ten=find(ten == 1);
    cont=zeros(2,4); % container for line i
    cont(1,1)=i;cont(2,1)=i; % id of line 
    cont(1,2)=str2double(C(indline(i)+ten)); % x start
    cont(1,3)=str2double(C(indline(i)+ten+2)); % y start
    cont(1,4)=str2double(C(indline(i)+ten+4)); % z start  
    cont(2,2)=str2double(C(indline(i)+ten+6)); % x end
    cont(2,3)=str2double(C(indline(i)+ten+8)); % y end
    cont(2,4)=str2double(C(indline(i)+ten+10)); % z end
    lines{i}=cont;
end
clear cont indline indnum ten linenum % delete garbage from workspace
%%%%%%%%%%%%%%%%%%%%%%%%%% end case line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% case polyline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get some markers and help variables
indpoly=strcmp('AcDbPolyline', C); % get line no. of LW polylines
polynum=sum(indpoly); % get total number of polylines
indpoly=find(indpoly == 1); % get line no. of polylines 
polylines=cell(polynum,1); % preallocate variable to increase speed
vertices=zeros(polynum,1); % detto
vertices=str2double(C(indpoly+2)); % the number of vertices is 2 lines after 'AcDbPolyline'
    for i=1:polynum % begin coordinate extraction for every single polyline, see !readme.m or dxf reference for details and group codes
        clear id xpoly ypoly zpoly  % clear to avoid error       
        null=strcmp('0',C(indpoly(i):(indpoly(i)+(4*vertices(i)+10)))); % find next 0 after last 10...=end of entity
        null=max(find(null == 1)); % max(null)=end of entity polyline(i)
        ten=strcmp('10',C((indpoly(i)+4):(indpoly(i)+null))); % find 10 in C (10 is group code for x-coords)
        ten=find(ten==1); % reshape ten
        NUM=str2double(C(indpoly(i):(indpoly(i)+null-1))); % get subset of numeric values of entity polyline(i), strings are nan's 
        xpoly=NUM(ten+5);ypoly=NUM(ten+7); % x- & y- coords
        threight=find(NUM(4:10)==38); % find '38' in NUM (38 is group code for z-coords) 
            if isempty(threight) % check out elevation
                zpoly=zeros(vertices(i),1); % elevation =0
            elseif threight~=0
                zpoly(1:vertices(i))=NUM(threight+4);zpoly=zpoly'; % get elevation if exists
            end          
        id(1:vertices(i))=i; % id of polyline
        polyline=[id' xpoly ypoly zpoly]; % create polylinesubset    
        seventy=find(NUM(3:6)==70); % find 70 (70 is group code for closed polylines)      
            if NUM(seventy+3)==1 % if polyline is closed...
                polyline(vertices(i)+1,:)=polyline(1,:); % ... add first row as last one
            end            
        polylines{i}=polyline; % save subset to array
        clear polyline % clear to avoid error
    end % end case polyline
clear i id indpoly null polynum seventy ten threight vertices 
clear xpoly ypoly zpoly NUM % delete garbage
%%%%%%%%%%%%%%%%%%%%%% end case polyline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% case 3d polyline %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get some markers and help variables
indpoly3d=strcmp('AcDb3dPolyline', C); % get line no. of 3d polylines
poly3dnum=sum(indpoly3d); % get total number of 3d polylines
indpoly3d=find(indpoly3d == 1); % detto 
vert3d=strcmp('AcDb3dPolylineVertex', C); % get line no. of vertices of 3d polylines
verttotal=sum(vert3d); % total amount of vertices
vert3d=find(vert3d == 1); % detto
poly3dlines=cell(poly3dnum,1); % preallocate 3dpolys

    for i=1:poly3dnum % begin coordinate extraction for every 3d polyline
        if i<poly3dnum
            idmax=max(find(vert3d<indpoly3d(i+1)));
        else  
            idmax=find(vert3d==max(vert3d));
        end            
            idmin=min(find(vert3d>indpoly3d(i)));
            sub=vert3d(idmin:1:idmax); % get indices of coords
            x3d=str2double(C(sub+2));
            y3d=str2double(C(sub+4));
            z3d=str2double(C(sub+6));
            id=zeros(1,length(x3d)); % to avoid error
            id(1:length(x3d))=i;
            ddpoly=zeros(length(sub),4);
            ddpoly=[id' x3d y3d z3d];            
            seventy=strcmp('70',C(indpoly3d+6:indpoly3d+11)); % find out if polyline is closed
            seventy=find(seventy==1);
            closed=str2double(C(indpoly3d(i)+seventy+6)); % get value of '70'+1 line
            closed=dec2bin(closed,8);closed=closed(:);closed=closed(8); % because group code '70' is binary coded
                        
            if closed=='1' % if polyline is closed 
                [x,y]=size(ddpoly);
                ddpoly(x+1,:)=ddpoly(1,:); % ... add first row as last one
            end
            
            poly3dlines{i}=ddpoly; % save subset to array
            clear ddpoly % clear to avoid error   


    end
clear closed ddpoly i id idmax idmin indpoly3d poly3dnum x3d y3d z3d
clear seventy sub vert3d verttotal x y z % delete garbage
%%%%%%%%%%%%%%%%%%%%%%%%%% end case 3d polyline %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% case 3d faces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get some markers and help variables
face3d=strcmp('AcDbFace', C); % get line no. of 3d faces
face3dnum=sum(face3d); % get total number of 3d faces
face3d=find(face3d == 1); % detto
faces3d=cell(face3dnum,1); % preallocate, +1 to enable resizing at the end

for i=1:face3dnum
    fac3d=zeros(4,4);
    fac3d(1,1)=i; % again(id,x-coord,y-coord,z_coords)
    fac3d(1,2)=str2double(C(face3d(i)+2));
    fac3d(1,3)=str2double(C(face3d(i)+4));
    fac3d(1,4)=str2double(C(face3d(i)+6));
    fac3d(2,1)=i;
    fac3d(2,2)=str2double(C(face3d(i)+8));
    fac3d(2,3)=str2double(C(face3d(i)+10));
    fac3d(2,4)=str2double(C(face3d(i)+12));
    fac3d(3,1)=i;
    fac3d(3,2)=str2double(C(face3d(i)+14));
    fac3d(3,3)=str2double(C(face3d(i)+16));
    fac3d(3,4)=str2double(C(face3d(i)+18));
    fac3d(4,1)=i;
    fac3d(4,2)=str2double(C(face3d(i)+20));
    fac3d(4,3)=str2double(C(face3d(i)+22));
    fac3d(4,4)=str2double(C(face3d(i)+24));
    
        if fac3d(4,2)==fac3d(3,2) % find out if 4th coord pair == 3rd
            fac3d=fac3d(1:3,:); % if so, delete 4th coord pair
        end
   faces3d{i}=fac3d; % save to cell array     
    
end
clear fac3d face3dnum face3d i % finished so far
%%%%%%%%%%%%%%%%%%%%%%%%%% end case 3d faces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% case circles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get some markers and help variables
cir=strcmp('AcDbCircle', C); % get line no. of 3d faces
cirnum=sum(cir); % get total number of 3d faces
cir=find(cir == 1); % detto
circles=cell(cirnum,1); % preallocate (id x y z radius)

for i=1:cirnum
    circ=zeros(1,5);
    threenine=strcmp('39',(C(cir(i):cir(i)+1))); % '39' is the group code of entity thickness
    threenine=find(threenine==1);
        if isempty (threenine) % get coordinates and radii
            circ=zeros(1,5);
            circ(1,1)=i;
            circ(1,2)=str2double(C(cir(i)+2)); % x-coord of center
            circ(1,3)=str2double(C(cir(i)+4)); % y-coord of center
            circ(1,4)=str2double(C(cir(i)+6)); % z-coord of center
            circ(1,5)=str2double(C(cir(i)+8)); % radius of circle
        else
            circ=zeros(1,5);
            circ(1,1)=i;
            circ(1,2)=str2double(C(cir(i)+4)); % x-coord of center
            circ(1,3)=str2double(C(cir(i)+6)); % y-coord of center
            circ(1,4)=str2double(C(cir(i)+8)); % z-coord of center
            circ(1,5)=str2double(C(cir(i)+10)); % radius of circle
        end

    circles{i}=circ; % save to array
    
end

clear threenine i cirnum circ cir  % delete garbage
%%%%%%%%%%%%%%%%%%%%%%%%%%% end case circles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rmpath(pathname); % remove path from the matlab search path











