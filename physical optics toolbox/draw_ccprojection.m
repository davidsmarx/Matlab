function hf = draw_ccprojection(haxes,cc,el)
% hf = draw_ccprojection(haxes,cc,edgelengths)

e0 = [0 0 0]'; e1 = cc(:,1); e2 = cc(:,2); e3 = cc(:,3);
% make face 1 = circular baseplate
NP = 21; % # number of points to make up circle
psilist = linspace(0,pi/2,NP);
R1z = make_rotmatrix(0,0,atan2(e1(2),e1(1)));
R1y = make_rotmatrix(0,acos(e1(3)),0);
face1 = e0;
for ii = 1:NP,
    tt = psilist(ii);
    rr = el(2)*el(3)*sqrt( (1+tan(tt).^2)/((el(2)*tan(tt)).^2 + el(3).^2) );
    face1 = [face1 ...
        rr*R1z*R1y*make_rotmatrix(0,0,psilist(ii))*R1y'*R1z'*e2];
end
%disp([face1(:,end) e3]);

% face 2
face2 = [e0 el(1)*e1 el(1)*e1+el(3)*e3 el(3)*e3];
% face 3
face3 = [e0 el(1)*e1 el(1)*e1+el(2)*e2 el(2)*e2];

if isempty(haxes),
    figure, haxes = axes;
else
    axes(haxes);
end
hf = fill(face1(1,:),face1(2,:),'b',...
    face2(1,:),face2(2,:),'r',...
    face3(1,:),face3(2,:),'g');
axis image, grid
xlabel('X'), ylabel('Y')

% h2 = figure;
% drawline = @(face,c) set(line(face(1,:),face(2,:),face(3,:)),'color',c);
% drawline([face1 face1(:,1)],'b')
% drawline([face2 face2(:,1)],'r')
% drawline([face3 face3(:,1)],'g')
% grid, axis image
% xlabel('X'), ylabel('Y'), zlabel('Z')

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3-d corner cube figure
figure, fill3(face1(1,:),face1(2,:),face1(3,:),'b',...
    face2(1,:),face2(2,:),face2(3,:),'r',...
    face3(1,:),face3(2,:),face3(3,:),'g'), grid
xlabel('X'), ylabel('Y'), zlabel('Z')

%          edge1theta = acos(e1(3));
%          edge1phi   = atan2(e1(2),e1(1))
%          DO ii = 1, NP
%             psi = (ii-1)*pi/2/NP
% 
%             cc%fM0(:,2+ii) = &
%                 MATMUL(Z_Rotation(edge1phi), &
%                  MATMUL(Y_Rotation(edge1theta), &
%                   MATMUL(Z_Rotation(psi), &
%                    MATMUL(Y_Rotation(-edge1theta), &
%                     MATMUL(Z_Rotation(-edge1phi), e2 ) ))))
%          ENDDO
% 
%          cc%fM0(:,2+NP+1) = e3
%          cc%fM0(:,2+NP+2) = e3 + e1
