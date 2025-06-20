function [nContact,nGap] = getThermalNetwork(sample,aX,aY,aZ,aR)
%用于生成热传导网络
%输入格式：[nContact,nGap]=f.run('fun/getThermalNetwork.m',d.GROUP.sample,d.mo.aX,d.mo.aY,d.mo.aZ,d.mo.aR);
m=size(sample,1);
X=aX(1:m);Y=aY(1:m);Z=aZ(1:m);R=aR(1:m);
Xmin=min(X)-min(R)/2;%以下四行用于确定网格边界
Xmax=max(X)+min(R)/2;
Ymin=min(Y)-min(R)/2;
Ymax=max(Y)+min(R)/2;
Zmin=min(Z)-min(R)/2;
Zmax=max(Z)+min(R)/2;
%dSide=0.4*mean(R);
gSide=4*max(R);%网格元胞边长
%gSide=6*max(R);
xNum=ceil((Xmax-Xmin)/gSide)+2;%以下两行划分网格
yNum=ceil((Ymax-Ymin)/gSide)+2;
zNum=ceil((Zmax-Zmin)/gSide)+2;
cellNum=xNum*yNum*zNum;%总的网格数
xIndex=floor((X-Xmin)/gSide)+2;%计算单元所在的网格
yIndex=floor((Y-Ymin)/gSide)+2;
zIndex=floor((Z-Zmin)/gSide)+2;
aIndex=xIndex+xNum*(yIndex-1+yNum*(zIndex-1));%将网格从1：ceelNum编号
[dcxIndex,dcyIndex,dczIndex]=ndgrid([-1,0,1],[-1,0,1],[-1,0,1]);%创建领域算子
cxIndex=xIndex'+dcxIndex(:);%记录中心网格的前后左右网格的行号
cyIndex=yIndex'+dcyIndex(:);%记录中心网格的前后左右网格的列号
czIndex=zIndex'+dczIndex(:);%记录中心网格的前后左右网格的列号
cIndex=cxIndex+xNum*(cyIndex-1+yNum*(czIndex-1));%单元-领域元胞矩阵
[sortaIndex,sortI]=sort(aIndex);%前为排序后的元胞列表，后为对应单元索引
[~,maxStepNum,~] = mode(aIndex);%单个元胞包含最多的单元数
cellBallList=zeros(maxStepNum,m);%初始化元胞-包含单元矩阵
%以下用于查找rowIndex
sortaIndexU=unique(sortaIndex);%sortaIndex向量重复元素的个数
rowIndex=zeros(length(sortaIndex),1);
for k=1:numel(sortaIndexU)
    ic=find(sortaIndex==sortaIndexU(k));%查找一组重复元素的位置
    ia=find(sortaIndex(ic));%给这组重复元素按顺序标号
    rowIndex(ic)=ia;
end
cellBallList(rowIndex+maxStepNum*(sortaIndex-1))=sortI;%获得元胞-包含单元矩阵
nBallGrid0=reshape(cellBallList(:,cIndex),27*maxStepNum,[])';
nBallGrid0=sort(nBallGrid0,2,'descend');
nBallGrid0(:,all(nBallGrid0==0,1))= [];%去掉全0列
zeroFilter=nBallGrid0==0;
selfFilter=nBallGrid0==(1:m)';
nBallGrid0(zeroFilter|selfFilter)=m+1;%将0和自身单元用虚单元填充
nBallGrid0=sort(nBallGrid0,2);
nBallGridAll=nBallGrid0(:,1:(end-1));%去除自身单元

%生成接触热传导矩阵
nBallGrid1=nBallGridAll;
Filter=nBallGrid1;%测量球-测量球过滤器
Filter(Filter==m+1)=0;
Filter(Filter~=0)=1; 
nBallGrid1(nBallGrid1==m+1)=m;
hij=((sqrt((X(nBallGrid1)-X(sample)).^2+(Y(nBallGrid1)-Y(sample)).^2+(Z(nBallGrid1)-Z(sample)).^2))-(R(nBallGrid1)+R(sample))).*Filter;
R2=zeros(size(nBallGrid1));
R2(hij<0)=1;
nBallGrid1=nBallGrid1.*R2;
nBallGrid1(nBallGrid1==0)=m+1;
nContact=sort(nBallGrid1,2,'ascend');
nContact(nContact==m+1)=0;
nContact(:,all(nContact==0,1))= [];%去掉全0列

%生成空气热传导矩阵
nBallGrid2=nBallGridAll;
nBallGrid2(nBallGrid2==m+1)=m;
R2=zeros(size(nBallGrid2));
Rij=0.5*(R(nBallGrid2).*R(sample))./(R(nBallGrid2)+R(sample)).*Filter;
R2(hij<Rij & hij>=0)=1;
nBallGrid2=nBallGrid2.*R2;
nBallGrid2(nBallGrid2==0)=m+1;
nGap=sort(nBallGrid2,2,'ascend');
nGap(nGap==m+1)=0;
nGap(:,all(nGap==0,1))= [];%去掉全0列


end