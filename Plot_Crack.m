clear
clc
count = 0;
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Choose appropriate value (Head) %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for AskTime = [2.0];
        %385 390 395 400 405 410 412 415 417 420 422 425 430 452]
        %239 261 274 284 294 311 322 331 340 347 358 374 387 401 410 419 437 457 472 490 502];
  
NumOfPlotEle = 294;
TimeTole=0.001; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hardening_strain=0.0001;
softening_strain=0.002;

Com_peak= -0.0075;
Com_crush =-0.3;

level_1=hardening_strain;
level_2=softening_strain;

level_3=10*softening_strain;
level_4=30*softening_strain;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hardening_strainUHPC=0.0001;
softening_strainUHPC=0.002;

Com_peakUHPC= -0.0075;
Com_crushUHPC =-0.3;

level_1UHPC=hardening_strain;
level_2UHPC=softening_strain;

level_3UHPC=10*softening_strain;
level_4UHPC=30*softening_strain;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UHPC = [1];
    
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Choose appropriate value (End) %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % node # and its coordinate
%NOdeCoor(totalnode#*4)=(node#,x,y,z)
NodeCoor=xlsread('Def_Coord.xlsx');
%%NodeCoor=xlsread('C:\z-PlottingCrack\NodeCoordinates.xlsx');
%total number of nodes = NumOfNode(1,1)
NumOfNode=size(NodeCoor(:,1));



% %element # and its nodes
%Ele(Total Ele # * 5)=(ele#,node1,node2,node3,node4)
Ele=xlsread('ele-try.xlsx');
% total element #
% total number of elements = NumOfEle(1,1)
NumOfEle=size(Ele(:,1));

figure()
title(AskTime)
hold on
% % for i=1:NumOfNode(1,1)
% % plot(NodeCoor(i,2), NodeCoor(i,3) ,'.','MarkerSize',0.1)
% % end

%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Read Principle direction (Head) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % node # and its coordinate
fid = fopen('hsv10');

    A=textscan(fid,'%n', 'commentStyle', '*');
    Ele_Num(1)=A{1};
    A=textscan(fid,'%n', 'commentStyle', 'endcurve');
    A=textscan(fid,'%n', 'commentStyle', '#');
    TimeSeri=textscan(fid,'%f %f', 'commentStyle', '*');
    Time_step=TimeSeri{1};
    Plane_Angle=TimeSeri{2};
    Ele_Num(1,2)=Plane_Angle(find(Time_step(:) > AskTime-TimeTole & Time_step(:) < AskTime+TimeTole, 1));

    for i=2: NumOfPlotEle
       
A=textscan(fid,'%n', 'commentStyle', 'endcurve')
count = count+1
Ele_Num(i,1)=A{1};

A=textscan(fid,'%n', 'commentStyle', 'endcurve');
A=textscan(fid,'%n', 'commentStyle', '#');
TimeSeri=textscan(fid,'%f %f', 'commentStyle', '*');
    Time_step=TimeSeri{1};
    Plane_Angle=TimeSeri{2};
    Ele_Num(i,2)=Plane_Angle(find(Time_step(:) > AskTime-TimeTole & Time_step(:) < AskTime+TimeTole, 1));
    end
    
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Read Principle direction (End) %%%%%%%
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Read hsv11 (Head) %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % node # and its coordinate
fid = fopen('hsv11');


    A=textscan(fid,'%n', 'commentStyle', '*');
%    Ele_Num(1)=A{1};
    A=textscan(fid,'%n', 'commentStyle', 'endcurve');
    A=textscan(fid,'%n', 'commentStyle', '#');
    TimeSeri=textscan(fid,'%f %f', 'commentStyle', '*');
    Time_step=TimeSeri{1};
    hsv11=TimeSeri{2};
    Ele_Num(1,3)=hsv11(find(Time_step(:) > AskTime-TimeTole & Time_step(:) < AskTime+TimeTole, 1))
    for i=2: NumOfPlotEle
A=textscan(fid,'%n', 'commentStyle', 'endcurve');
%Ele_Num(i,1)=A{1};
A=textscan(fid,'%n', 'commentStyle', 'endcurve');
A=textscan(fid,'%n', 'commentStyle', '#');
TimeSeri=textscan(fid,'%f %f', 'commentStyle', '*');
    Time_step=TimeSeri{1};
    hsv11=TimeSeri{2};
    Ele_Num(i,3)=hsv11(find(Time_step(:) > AskTime-TimeTole & Time_step(:) < AskTime+TimeTole, 1));
    end

    
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Read hsv11 (End) %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Read hsv12 (Head) %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % node # and its coordinate
fid = fopen('hsv12');

    A=textscan(fid,'%n', 'commentStyle', '*');
%    Ele_Num(1)=A{1};
    A=textscan(fid,'%n', 'commentStyle', 'endcurve');
    A=textscan(fid,'%n', 'commentStyle', '#');
    TimeSeri=textscan(fid,'%f %f', 'commentStyle', '*');
    Time_step=TimeSeri{1};
    hsv12=TimeSeri{2};
    Ele_Num(1,4)=hsv12(find(Time_step(:) > AskTime-TimeTole & Time_step(:) < AskTime+TimeTole, 1));

    for i=2: NumOfPlotEle
A=textscan(fid,'%n', 'commentStyle', 'endcurve');
%Ele_Num(i,1)=A{1};
A=textscan(fid,'%n', 'commentStyle', 'endcurve');
A=textscan(fid,'%n', 'commentStyle', '#');
TimeSeri=textscan(fid,'%f %f', 'commentStyle', '*');
    Time_step=TimeSeri{1};
    hsv12=TimeSeri{2};
    Ele_Num(i,4)=hsv12(find(Time_step(:) > AskTime-TimeTole & Time_step(:) < AskTime+TimeTole, 1));
    end
    
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Read hsv12 (End) %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Plot Crack  (Head) %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:NumOfPlotEle

    
    Ele_Num(i,1)
    entryy = ismember(Ele_Num(i,1),UHPC);
    

if  entryy == 0;
    
    % find position of the corresponding element number: find( Ele(:,1) == Ele_Num(i) )
    CorrespondEN = find( Ele(:,1) == Ele_Num(i,1) );
    Node1= Ele(CorrespondEN,2);
    Node2= Ele(CorrespondEN,3);
    Node3= Ele(CorrespondEN,4);
    Node4= Ele(CorrespondEN,5);
    Node1_XY(1:2)=NodeCoor(find(NodeCoor(:,1) == Node1),2:3);
    Node2_XY(1:2)=NodeCoor(find(NodeCoor(:,1) == Node2),2:3);
    Node3_XY(1:2)=NodeCoor(find(NodeCoor(:,1) == Node3),2:3);
    Node4_XY(1:2)=NodeCoor(find(NodeCoor(:,1) == Node4),2:3);
    
    %coordinate of the center node in an element
    
    Coord_Cen(1:2)=(Node1_XY+Node2_XY+Node3_XY+Node4_XY)/4;
    
    %average side length
    
    xxx=Node1_XY-Node2_XY;
    Len(1)=(xxx(1)^2+xxx(2)^2)^0.5;
    xxx=Node2_XY-Node3_XY;
    Len(2)=(xxx(1)^2+xxx(2)^2)^0.5;
    xxx=Node3_XY-Node4_XY;
    Len(3)=(xxx(1)^2+xxx(2)^2)^0.5;
    xxx=Node4_XY-Node1_XY;
    Len(4)=(xxx(1)^2+xxx(2)^2)^0.5;
    Ave_Len=sum(Len)/4;
    
    % hsv11
    
    if ((Ele_Num(i,3) >= hardening_strain) || (Ele_Num(i,3) <= Com_peak));
    Coor_End_1=[Coord_Cen(1)-Ave_Len/2*sin(Ele_Num(i,2)), Coord_Cen(2)+cos(Ele_Num(i,2))*Ave_Len/2];
    Coor_End_2=[Coord_Cen(1)+Ave_Len/2*sin(Ele_Num(i,2)), Coord_Cen(2)-cos(Ele_Num(i,2))*Ave_Len/2];
    
    X_Coor=[Coor_End_1(1,1) Coor_End_2(1,1)];
    Y_Coor=[Coor_End_1(1,2) Coor_End_2(1,2)];
    
%   figure(1)

    % Crack
    
    if ((Ele_Num(i,3) >= hardening_strain) && (Ele_Num(i,3) < level_2));
    line(X_Coor,Y_Coor,'Color','r','LineStyle',':','LineWidth',2);
    z=11;
    elseif ((Ele_Num(i,3) >= level_2) && (Ele_Num(i,3) < level_3));
    line(X_Coor,Y_Coor,'Color','r','LineStyle','--','LineWidth',2);
    z=12;
    elseif ((Ele_Num(i,3) >= level_3) && (Ele_Num(i,3) < level_4));
    line(X_Coor,Y_Coor,'Color','r','LineWidth',4);
    z=13;
    elseif ((Ele_Num(i,3) >= level_4) );
    line(X_Coor,Y_Coor,'Color','r','LineWidth',8);
    z=14;
    
    % Crush
    
    elseif ((Ele_Num(i,3) <= Com_peak) && (Ele_Num(i,3) > Com_crush));
    line(X_Coor,Y_Coor,'Color','k','LineWidth',3);
    plot(X_Coor,Y_Coor,'.','markersize',15);
    z=15;
    elseif ((Ele_Num(i,3) <= Com_crush) );
    line(X_Coor,Y_Coor,'Color','k','LineWidth',4);
    plot(X_Coor,Y_Coor,'.','markersize',40);
    z=16;
    end
    
    end
    
    
    
    % hsv12
    
    if ((Ele_Num(i,4) >= hardening_strain) || (Ele_Num(i,4) <= Com_peak));
    Coor_End_1_hsv12=[Coord_Cen(1)+Ave_Len/2*cos(Ele_Num(i,2)), Coord_Cen(2)+sin(Ele_Num(i,2))*Ave_Len/2];
    Coor_End_2_hsv12=[Coord_Cen(1)-Ave_Len/2*cos(Ele_Num(i,2)), Coord_Cen(2)-sin(Ele_Num(i,2))*Ave_Len/2];
    
    X_Coor_hsv12=[Coor_End_1_hsv12(1,1) Coor_End_2_hsv12(1,1)];
    Y_Coor_hsv12=[Coor_End_1_hsv12(1,2) Coor_End_2_hsv12(1,2)];
    
    % figure(1)
    % Crack
    
    if ((Ele_Num(i,4) >= hardening_strain) && (Ele_Num(i,4) < level_2));
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','m','LineStyle',':','LineWidth',2);
    z=21;
    elseif ((Ele_Num(i,4) >= level_2) && (Ele_Num(i,4) < level_3));
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','m','LineStyle','--','LineWidth',2);
    z=22;
    elseif ((Ele_Num(i,4) >= level_3) && (Ele_Num(i,4) < level_4));
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','m','LineWidth',4);
    z=23;
    elseif ((Ele_Num(i,4) >= level_4) );
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','m','LineWidth',8); 
    z=24;
    %Crush
    elseif ((Ele_Num(i,4) <= Com_peak) && (Ele_Num(i,4) > Com_crush));
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','g','LineWidth',5);
    plot(X_Coor_hsv12,Y_Coor_hsv12,'.','markersize',15);
    z=25;
    elseif ((Ele_Num(i,3) <= Com_crush) );
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','g','LineWidth',5);
    plot(X_Coor_hsv12,Y_Coor_hsv12,'.','markersize',15);
    z=26;
    end
    
    end
    
else 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % find position of the corresponding element number: find( Ele(:,1) == Ele_Num(i) )
    CorrespondEN = find( Ele(:,1) == Ele_Num(i,1) );
    Node1= Ele(CorrespondEN,2);
    Node2= Ele(CorrespondEN,3);
    Node3= Ele(CorrespondEN,4);
    Node4= Ele(CorrespondEN,5);
    Node1_XY(1:2)=NodeCoor(find(NodeCoor(:,1) == Node1),2:3);
    Node2_XY(1:2)=NodeCoor(find(NodeCoor(:,1) == Node2),2:3);
    Node3_XY(1:2)=NodeCoor(find(NodeCoor(:,1) == Node3),2:3);
    Node4_XY(1:2)=NodeCoor(find(NodeCoor(:,1) == Node4),2:3);
    
    %coordinate of the center node in an element
    
    Coord_Cen(1:2)=(Node1_XY+Node2_XY+Node3_XY+Node4_XY)/4;
    
    %average side length
    
    xxx=Node1_XY-Node2_XY;
    Len(1)=(xxx(1)^2+xxx(2)^2)^0.5;
    xxx=Node2_XY-Node3_XY;
    Len(2)=(xxx(1)^2+xxx(2)^2)^0.5;
    xxx=Node3_XY-Node4_XY;
    Len(3)=(xxx(1)^2+xxx(2)^2)^0.5;
    xxx=Node4_XY-Node1_XY;
    Len(4)=(xxx(1)^2+xxx(2)^2)^0.5;
    Ave_Len=sum(Len)/4;
    
    % hsv11
    
    if ((Ele_Num(i,3) >= hardening_strainUHPC) || (Ele_Num(i,3) <= Com_peakUHPC));
    Coor_End_1=[Coord_Cen(1)-Ave_Len/2*sin(Ele_Num(i,2)), Coord_Cen(2)+cos(Ele_Num(i,2))*Ave_Len/2];
    Coor_End_2=[Coord_Cen(1)+Ave_Len/2*sin(Ele_Num(i,2)), Coord_Cen(2)-cos(Ele_Num(i,2))*Ave_Len/2];
    
    X_Coor=[Coor_End_1(1,1) Coor_End_2(1,1)];
    Y_Coor=[Coor_End_1(1,2) Coor_End_2(1,2)];
    
%   figure(1)

    % Crack
    
    if ((Ele_Num(i,3) >= hardening_strainUHPC) && (Ele_Num(i,3) < level_2UHPC));
    line(X_Coor,Y_Coor,'Color','g','LineStyle',':','LineWidth',2);
    elseif ((Ele_Num(i,3) >= level_2UHPC) && (Ele_Num(i,3) < level_3UHPC));
    line(X_Coor,Y_Coor,'Color','g','LineStyle','--','LineWidth',2);
    elseif ((Ele_Num(i,3) >= level_3UHPC) && (Ele_Num(i,3) < level_4UHPC));
    line(X_Coor,Y_Coor,'Color','g','LineWidth',4);
    elseif ((Ele_Num(i,3) >= level_4UHPC) );
    line(X_Coor,Y_Coor,'Color','g','LineWidth',8);
    
    % Crush
    
    elseif ((Ele_Num(i,3) <= Com_peakUHPC) && (Ele_Num(i,3) > Com_crushUHPC));
    line(X_Coor,Y_Coor,'Color','m','LineWidth',3);
    plot(X_Coor,Y_Coor,'.','markersize',15);
    elseif ((Ele_Num(i,3) <= Com_crushUHPC) );
    line(X_Coor,Y_Coor,'Color','m','LineWidth',4);
    plot(X_Coor,Y_Coor,'.','markersize',40);
    end
    
    end
    
    
    
    % hsv12
    
    if ((Ele_Num(i,4) >= hardening_strainUHPC) || (Ele_Num(i,4) <= Com_peakUHPC));
    Coor_End_1_hsv12=[Coord_Cen(1)+Ave_Len/2*cos(Ele_Num(i,2)), Coord_Cen(2)+sin(Ele_Num(i,2))*Ave_Len/2];
    Coor_End_2_hsv12=[Coord_Cen(1)-Ave_Len/2*cos(Ele_Num(i,2)), Coord_Cen(2)-sin(Ele_Num(i,2))*Ave_Len/2];
    
    X_Coor_hsv12=[Coor_End_1_hsv12(1,1) Coor_End_2_hsv12(1,1)];
    Y_Coor_hsv12=[Coor_End_1_hsv12(1,2) Coor_End_2_hsv12(1,2)];
    
    % figure(1)
    % Crack
    
    if ((Ele_Num(i,4) >= hardening_strainUHPC) && (Ele_Num(i,4) < level_2UHPC));
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','g','LineStyle',':','LineWidth',2);
    elseif ((Ele_Num(i,4) >= level_2UHPC) && (Ele_Num(i,4) < level_3UHPC));
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','g','LineStyle','--','LineWidth',2);
    elseif ((Ele_Num(i,4) >= level_3UHPC) && (Ele_Num(i,4) < level_4UHPC));
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','g','LineWidth',4);
    elseif ((Ele_Num(i,4) >= level_4UHPC) );
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','g','LineWidth',8); 
    %Crush
    elseif ((Ele_Num(i,4) <= Com_peakUHPC) && (Ele_Num(i,4) > Com_crushUHPC));
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','m','LineWidth',5);
    plot(X_Coor_hsv12,Y_Coor_hsv12,'.','markersize',15);
    elseif ((Ele_Num(i,3) <= Com_crushUHPC) );
    line(X_Coor_hsv12,Y_Coor_hsv12,'Color','m','LineWidth',5);
    plot(X_Coor_hsv12,Y_Coor_hsv12,'.','markersize',15);
    end
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
end
end


end

%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Plot Crack  (End) %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 plot(NodeCoor(:,2),NodeCoor(:,3),'bx')
 title('Cracking of Regular Concrete Box Girder with Regular Concrete Joint')
 set(gca,'xtick',[], 'xticklabel',{})
 set(gca,'ytick',[], 'yticklabel',{})

% 
% 
% %%%%%%%%%%%%%%%%%%% Plot all the elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 figure(1)
 hold on
 for j=1:NumOfEle(1,1)
     for i=1:3
         X1=[NodeCoor(find(NodeCoor(:,1) == Ele(j,i+1)),2), NodeCoor(find(NodeCoor(:,1) == Ele(j,i+2)),2)];
         Y1=[NodeCoor(find(NodeCoor(:,1) == Ele(j,i+1)),3), NodeCoor(find(NodeCoor(:,1) == Ele(j,i+2)),3)];
         figure(1)
         line(X1,Y1,'Color','b','LineWidth',1)
     end
 
     X2=[NodeCoor(find(NodeCoor(:,1) == Ele(j,5)),2), NodeCoor(find(NodeCoor(:,1) == Ele(j,2)),2)];
     Y2=[NodeCoor(find(NodeCoor(:,1) == Ele(j,5)),3), NodeCoor(find(NodeCoor(:,1) == Ele(j,2)),3)];
     figure(1)
     line(X2,Y2,'Color','b','LineWidth',1)
 end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% 
% saveas(figure(1),'Cracks.jpg','jpg');