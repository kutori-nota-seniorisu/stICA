function twitchCurves = twitchCurve(twitchFrames,actArea)

twitchCurves = cell(1,size(actArea,2));

for na = 1:size(actArea,2)
    tmpCurves_all = [];
    for np = 1:size(actArea{na},1)
        % 抽出一个激活像素点处的数值
        tmp = twitchFrames(actArea{na}(np,1),actArea{na}(np,2),:);
        % 降维
        tmp = reshape(tmp,1,[]);
        tmpCurves_all = [tmpCurves_all;tmp];
    end
    twitchCurves{na} = mean(tmpCurves_all, 1);
end


