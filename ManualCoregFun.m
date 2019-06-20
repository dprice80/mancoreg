function [] = ManualCoregFun(Filename)
% Manual Fiducial Coregistration
% Darren Price - 25/06/14
% Inputs
% Filename = filename of resliced NII file with 1mm isotropic resolution.
% Input file should ideally be recogistered to MNI space using a rigid body
% transformation and resliced (reslicing is very important).

if nargin < 2
    S.blurval = 1.2;
    S.overwrite = 1;
    S.Ang = 0; % (0->360) Lighting angle in degrees.
    S.ShaderFade = 0.2; % (0->1) Lower = smoother, less defined; Higher = more defined, but voxelated
    S.DepthFade = 0.6; % (0->1) Lower = less depth fade
    S.save = 1;
end

blurval = S.blurval;
overwrite = S.overwrite;
Ang = (S.Ang+90)*360/(2*pi);
ShaderFade = S.ShaderFade;
DepthFade = S.DepthFade;

if ~isfield(S,'save')
    S.save = 1;
end

[savedir,fname,~] = fileparts(Filename);

if exist([savedir '/' fname '_fid.mat'],'file') || exist([savedir '/Flagged.mat'],'file')
    disp(['Subject Completed or Flagged: ' Filename])
    if overwrite == 0
        disp(['Subject in progress: ' Filename])
        return
    end
end

NII = load_nii(Filename);

origin = NII.hdr.hist.originator(1:3); % pull out origin
scale = NII.hdr.dime.pixdim; % pull out dimensions
dims = NII.hdr.dime.dim; %#ok<NASGU>
hdr = NII.hdr; %#ok<NASGU> % save the header

T = [NII.hdr.hist.srow_x; NII.hdr.hist.srow_y; NII.hdr.hist.srow_z; [0 0 0 1]];

niiOrig = NII.img;
% H = vol3d('CData',niiOrig);
correct = 0;
threshlim = 4;


while correct == 0;
    close all
    
    nii = double(niiOrig);
    [~, bins] = hist(nii(nii>0),50);
    thresh = bins(threshlim);
    nii(nii < thresh) = 2;
    niiL = zeros(size(nii,2),size(nii,3));
    niiR = zeros(size(nii,2),size(nii,3));
    niiN = zeros(size(nii,1),size(nii,3));
    
    % Left
    for yi = 1:size(nii,2)
        for zi = 1:size(nii,3)
            findL = find(nii(:,yi,zi)>thresh,1,'first');
            if ~isempty(findL) % if its empty then skip row entirely
                niiL(yi,zi)  = findL;
            else
                niiL(yi,zi) = 0;
            end
        end
    end
    
    % Right
    for yi = 1:size(nii,2)
        for zi = 1:size(nii,3)
            findR = find(nii(:,yi,zi)>thresh,1,'last');
            if ~isempty(findR) % if its empty then skip row entirely
                niiR(yi,zi) = findR;%nii(findR,yi,zi); % nii(findR,yi,zi); % this is the skin surface
            else
                niiR(yi,zi) = 0;
            end
        end
    end
    
    % Nasion
    for xi = 1:size(nii,1)
        for zi = 1:size(nii,3)
            findN = find(nii(xi,:,zi)>thresh,1,'last');
            if ~isempty(findN) % if its empty then skip row entirely
                niiN(xi,zi) = findN;%nii(findN,yi,zi); % nii(findN,yi,zi); % this is the skin surface
            else
                niiN(xi,zi) = 0;
            end
        end
    end
    %%
    cmap1 = colormap('gray');
    monpos = get(0,'MonitorPositions');
    % figsize = [monpos(3)/4 monpos(4)/4 monpos(3)/2 monpos(4)/2];
    
    % Plot the data
    
    LAng = -[sin(Ang);cos(Ang)];
    niiL(niiL==0) = max(niiL(:));
    
%     figure(100)
%     imshow('/imaging/camcan/sandbox/projects/megcoregister/Scripts/rpa.JPG','InitialMagnification',25)
%     set(gcf,'position',[200 monpos(3)/3-200 200 400])
    [yL, zL] = CalcGrads(niiL,LAng,cmap1,blurval,ShaderFade,DepthFade);
%     figure(100)
%     imshow('/imaging/camcan/sandbox/projects/megcoregister/Scripts/rpa.JPG','InitialMagnification',25)
%     set(gcf,'position',[200 monpos(3)/3-200 200 400])
    [yR, zR] = CalcGrads(max(niiR(:))-niiR,LAng,cmap1,blurval,ShaderFade,DepthFade);
%     figure(100)
%     imshow('/imaging/camcan/sandbox/projects/megcoregister/Scripts/Nose.JPG','InitialMagnification',25)
%     set(gcf,'position',[200 monpos(3)/3-200 200 400])
    [xN, zN] = CalcGrads(max(niiN(:))-niiN,LAng,cmap1,blurval,ShaderFade,DepthFade);
    
    
    yL = round(yL(end)); zL = round(zL(end)); xL = round(niiL(yL,zL));
    yR = round(yR(end)); zR = round(zR(end)); xR = round(niiR(yR,zR));
    xN = round(xN(end)); zN = round(zN(end)); yN = round(niiN(xN,zN));
    
    M = [xL yL zL;xR yR zR;xN yN zN];
    
    % Convert to mm.
    Mmm = M-repmat(origin,3,1);
    
    % X Direction
    % X L
    MarkerSize = 20;
    figure(1); set(gcf,'position',monpos)
    subplot(3,3,1)
    imagesc(squeeze(nii(xL,:,:))'); axis equal; hold on
    set(gca,'XDir','normal','YDir','normal'); colormap('gray')
    plot(yL,zL,'+w','MarkerSize',MarkerSize); title('Left')
    % X R
    subplot(3,3,2)
    imagesc(squeeze(nii(xR,:,:))'); axis equal; hold on
    set(gca,'XDir','normal','YDir','normal'); colormap('gray')
    plot(yR,zR,'+w','MarkerSize',MarkerSize); title('Right')
    % X N
    subplot(3,3,3)
    imagesc(squeeze(nii(xN,:,:))'); axis equal; hold on
    set(gca,'XDir','normal','YDir','normal'); colormap('gray')
    plot(yN,zN,'+w','MarkerSize',MarkerSize); title('Nasion')
    
    % Y Direction
    % Y L
    subplot(3,3,4)
    imagesc(squeeze(nii(:,yL,:))'); axis equal; hold on
    set(gca,'XDir','normal','YDir','normal'); colormap('gray')
    plot(xL,zL,'+w','MarkerSize',MarkerSize);
    title(['x=' num2str(Mmm(1,1)) ' y=' num2str(Mmm(1,2)) ' z=' num2str(Mmm(1,3)) 'mm']) % put location in mm on plots
    % Y R
    subplot(3,3,5)
    imagesc(squeeze(nii(:,yR,:))'); axis equal; hold on
    set(gca,'XDir','normal','YDir','normal'); colormap('gray')
    plot(xR,zR,'+w','MarkerSize',MarkerSize);
    title(['x=' num2str(Mmm(2,1)) ' y=' num2str(Mmm(2,2)) ' z=' num2str(Mmm(2,3)) 'mm'])
    % Y N
    subplot(3,3,6)
    imagesc(squeeze(nii(:,yN,:))'); axis equal; hold on
    set(gca,'XDir','normal','YDir','normal'); colormap('gray')
    plot(xN,zN,'+w','MarkerSize',MarkerSize);
    title(['x=' num2str(Mmm(3,1)) ' y=' num2str(Mmm(3,2)) ' z=' num2str(Mmm(3,3)) 'mm'])
    
    % Z Direction
    % Z L
    subplot(3,3,7)
    imagesc(squeeze(nii(:,:,zL))'); axis equal; hold on
    set(gca,'XDir','normal','YDir','normal'); colormap('gray')
    plot(xL,yL,'+w','MarkerSize',MarkerSize);
    % Z R
    subplot(3,3,8)
    imagesc(squeeze(nii(:,:,zR))'); axis equal; hold on
    set(gca,'XDir','normal','YDir','normal'); colormap('gray')
    plot(xR,yR,'+w','MarkerSize',MarkerSize);
    % Z N
    subplot(3,3,9)
    imagesc(squeeze(nii(:,:,zN))'); axis equal; hold on
    set(gca,'XDir','normal','YDir','normal'); colormap('gray')
    plot(xN,yN,'+w','MarkerSize',MarkerSize);
    
    pbL = uicontrol('Style','pushbutton','string','Update Left','position',[20 20 80 30],'Callback',@UpdateCoordsCallbackL);
    pbR = uicontrol('Style','pushbutton','string','Update Right','position',[100 20 80 30],'Callback',@UpdateCoordsCallbackR);
    pbN = uicontrol('Style','pushbutton','string','Update Nasion','position',[180 20 80 30],'Callback',@UpdateCoordsCallbackN);
    
    pause(0.1)
    
    %% Decide what to do next save/rerun/flag
    choice = menu('Happy?','YES GIMME MORE!','NO, Try Again','NO Flag This Subject','Save and QUIT');
    
    % update coords if user updated
    udL = get(pbL,'UserData');
    if ~isempty(udL)
        M(1,1) = round(udL.m);
        M(1,2) = round(udL.n);
    end
    
    udR = get(pbR,'UserData');
    if ~isempty(udR)
        M(2,1) = round(udR.m);
        M(2,2) = round(udR.n);
    end
    
    udN = get(pbN,'UserData');
    if ~isempty(udN)
        M(3,2) = round(udN.m);
        M(3,3) = round(udN.n);
    end
    
    % Convert to mm.
    Mmm = M-repmat(origin,3,1);

    
    switch choice
        case 1 % save the data and clear up the mess
            correct = 1;
            save([savedir '/%s' fname '_fid.mat'],'Mmm')
        case 2 % do nothing to rerun while-loop
            
        case 3 % flag up this subject
            correct = 3;
            save([savedir '/Flagged.mat'],'')
        case 4
            correct = 1;
            save([savedir '/' fname '_fid.mat'],'Mmm')
            Quitting = 1; %#ok<NASGU>
    end
    
    if any(M(:) == 0)
        correct = 0;
        disp('Zeros detected: Retrying subject.')
    end
    close all
end

return

    function [m, n] = CalcGrads(niiIn,L,cmap,blurval,ShaderFade,DepthFade)
        [gx, gy] = gradient(niiIn); % get gradient of the surface for shading
        gxr = reshape(gx,1,numel(gx)); % turn in to a vector so we can take the dot prod with the lighting direction vector L
        gyr = reshape(gy,1,numel(gy));
        B = reshape([gxr; gyr]'*L,size(gx)); % dot prod (find angle between surface and light) and return reshaped matrix
        if blurval > 1 % gaussian blurring, get rid of voxelation (not really necessary)
            h = fspecial('gaussian', size(gx), blurval); % create a gaussian blurring kernel
            B = conv2(B,h,'same'); % convolve that kernel with the gradient image.
        end
        B = sin(B); % find sin of the angle of incidence
        niiIn = niiIn-min(niiIn(:)); % rescale input matrix
        niiIn = niiIn./max(max(niiIn(:,1:70)));
        niiIn(niiIn > 1) = 1;
        % Create depth mask
        niiInShad = max(niiIn(:))-niiIn; % invert input matrix
        niiInShad(niiInShad < DepthFade) = DepthFade; % set depth limit
        niiInShad = niiInShad-min(niiInShad(:)); % reset 0
        niiInShad = niiInShad./max(niiInShad(:));  % rescale
        BS = B.*niiInShad;
        BS = BS./max(BS(:)); % rescale gradient field
        niiShaded = niiInShad+BS.*ShaderFade; % apply depth shading
        clim = max(max(niiShaded(:,1:70)));
        figure(1);
        monpos = get(0,'MonitorPositions');
        set(gcf,'position',[400 monpos(3)/3-200 600 600])
        imagesc(niiShaded'); axis equal;
        set(gca,'XDir','normal','YDir','normal')
        % cmap = (cmap.^2)./(max(max(cmap.^2)))
        colormap(cmap);
        caxis([0 clim+0.1])
        title('LEFT or SPACEBAR = select, RIGHT = zoom, DOUBLE-CLICK RIGHT = zoom out')
        [m,n] = ginput2(1,false,'.r');
        hold on; plot(m(end),n(end),'+','MarkerSize',30); hold off; pause(1);
    end

    function UpdateCoordsCallbackL(hObject, eventdata, handles)
        % hObject    handle to pushbutton1 (see GCBO)
        % eventdata  reserved - to be defined in a future version of MATLAB
        % handles    structure with handles and user data (see GUIDATA)
        figure(1);subplot(3,3,7)
        [m,n] = ginput2(1,false,'.r');
        hold on; plot(m(end),n(end),'+','MarkerSize',30); hold off; pause(1);
        Out.m = m;
        Out.n = n;
        disp('Coordinates Updated')
        set(hObject,'UserData',Out)
    end

    function UpdateCoordsCallbackR(hObject, eventdata, handles)
        figure(1);subplot(3,3,7)
        [m,n] = ginput2(1,false,'.r');
        hold on; plot(m(end),n(end),'+','MarkerSize',30); hold off; pause(1);
        Out.m = m;
        Out.n = n;
        disp('Coordinates Updated')
        set(hObject,'UserData',Out)
    end

    function UpdateCoordsCallbackN(hObject, eventdata, handles)
        figure(1);subplot(3,3,7)
        [m,n] = ginput2(1,false,'.r');
        hold on; plot(m(end),n(end),'+','MarkerSize',30); hold off; pause(1);
        Out.m = m;
        Out.n = n;
        disp('Coordinates Updated')
        set(hObject,'UserData',Out)
    end

end