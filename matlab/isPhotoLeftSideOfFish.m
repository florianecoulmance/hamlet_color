function leftFlag = isPhotoLeftSideOfFish(short_name)


% Right or left side?
dashLoc = strfind(short_name,'-');
leftFlag = strcmp(short_name(dashLoc(1)+1),'l');