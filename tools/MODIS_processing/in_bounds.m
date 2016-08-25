function blIn = in_bounds(UL, LR, crd)

%crd = [W, E, S, N];
%UL = [x,y] for upper-left corner
%LR = [x,y] for lower-right corner

blIn = 1;

if LR(1) < crd(1) %tile is too far West
	blIn = 0; 
elseif UL(1) > crd(2) %tile is too far East
	blIn = 0;     
elseif UL(2) < crd(3) %tile is too far South
	blIn = 0;     
elseif LR(2) > crd(4) %tile is too far North
	blIn = 0;     
end


% if (LR(1) < crd(1) && UL(1) < crd(1)) %tile is too far West
% 	blIn = 0; 
% elseif (LR(1) > crd(2) && UL(1) > crd(2)) %tile is too far East
% 	blIn = 0;     
% elseif (LR(2) < crd(3) && UL(2) < crd(3)) %tile is too far South
% 	blIn = 0;     
% elseif (LR(2) > crd(4) && UL(2) > crd(4)) %tile is too far North
% 	blIn = 0;     
% end