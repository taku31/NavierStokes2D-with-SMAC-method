% Solving 2D navier stokes equation with SMAC method
% Copyright (C) 2019  T.Nakabayashi
% Released under the MIT license http://opensource.org/licenses/mit-license.php

% 初期化
clear all;

% グローバル変数宣言
global dt ddt nx ny dx dy ddx ddx2 ddy ddy2 re

% パラメーター
n = 25;% 格子数
nx = 6 * n;% ｘ方向格子数
ny = 2 * n;% % ｙ方向格子数
loop = 20000;% ステップ数
re = 100;% レイノルズ数
dt = 0.02;% タイムステップ

% 配列の確保
p = zeros(nx + 2, ny + 2);
u = zeros(nx + 1, ny + 2);
v = zeros(nx + 2, ny + 1);
phi = zeros(nx + 2, ny + 2);% 補正圧力
up = zeros(nx + 1,ny + 2);% 予測速度
vp = zeros(nx + 2,ny + 1);% 予測速度
divup = zeros(nx + 2, ny + 2);% 予測速度の発散
divu = zeros(nx + 2, ny + 2);% 連続の式チェック用
psi = zeros(nx + 1, ny + 2);% 流れ関数
uu = zeros(nx+2,ny+2);
vv = zeros(nx+2,ny+2);

% 除算数の削減
dx = 5 / n;% 格子幅
dy = dx;
ddx = 1 / dx;
ddy = ddx;
ddx2 = ddx * ddx;
ddy2 = ddy * ddy;
ddt = 1 / dt;

% タイムステップの大きさ確認
dt = min(dt, 0.25 * dx);% 1時間のステップで流体が移流によって飛び出さないようにする
dt = min(dt, 0.2 * re * dx * dx);% 拡散の影響の考慮

% 初期条件の代入
u(:,ny+2) = 1;% 全領域を１にすることで、初タイムステップでできる限り連続の方程式を満たすようする

% 境界条件の設定
un = 1;
uw = 1;% 流入口 
us = 1;
ue = 0;% 流出口
vn = 0;
vw = 0;% 流入口
vs = 0;
ve = 0;% 流出口

% 障害物位置の定義
object = zeros(nx + 2, ny + 2);% 圧力格子ベースで障害物を定義する。
center = [(nx + 2) / 6, (ny + 3) / 2];
object = DefineObjectArea(object, center);

for ita = 1 : loop
    
    disp(ita);
    
    % 修正流速・修正圧力u,v,pへ境界条件を適用
    u = BoundaryConditionU(u, ue, uw, us, un);
    v = BoundaryConditionV(v, ve, vw, vs, vn);
    p = BoundaryConditionP(p);
    
    % 仮流速upの計算
    up = ProvisionalVelocityU(u, v, p, up);
    
    % 仮流速upへ境界条件を適用
    up = BoundaryConditionU(up, ue, uw, us, un);
    
    % 仮流速vpの計算
    vp = ProvisionalVelocityV(u, v, p, vp);
    
    % 仮流速vpへ境界条件を適用
    vp = BoundaryConditionV(vp, ve, vw, vs, vn);
    
    % 仮流速が連続の方程式を満たしているか確認
    [div, divup] = CheackContinuityEquation(up, vp, divup);
    
    % 圧力のポアソン方程式を解く
    eps = 10^(- 8);
    maxitr = nx * ny * 2;% 反復回数。収束させるためにはこのぐらい必要。
    alpha = 1.7;% 緩和係数
    phi = PoissonSolver(alpha, phi, eps, maxitr, divup, nx, ny, ddt, ddx2, ddy2);
    
    % 仮速度・圧力の修正
    [u, v, p] = ModifyVP(up, vp, u, v, p, phi);
    
    % 修正流速が連続の式の満足度チェック
    [div, divu] = CheackContinuityEquation(u, v, divu);
    
    % 障害物の配置(セル内で占有する割合によって値を決定する。)
    for i = 1 : nx + 2
        for j = 1 : ny + 2
            if object(i, j) == 1% 障害物内部ならば値をゼロに置く
                u(i - 1, j) = 0;
                v(i, j - 1) = 0;
            elseif object(i, j) == 2% 障害物境界ならば値の半分
                u(i - 1, j) = 0.5 * u(i - 1, j);
                v(i, j - 1) = 0.5 * v(i, j - 1);
            end
        end
    end
    
    % 圧力格子位置での速度を求める。
    [uu, vv] = VelocityInterpolate(u, v, uu, vv);
    
    % 結果の描画
    vis_contour('u.gif', ita, uu, 0, 1.5, 1)
    %vis_contour('v.gif', ita, vv, -0.6, 0.6, 2)
    %vis_vector('vec.gif', ita, uu, vv, 3)
    
end

function[] = vis_contour(filename, timestep, u, maxrange, minrange, fignum)
% スカラー場の可視化%
% Input
% ----------
% filename : text
%   出力gifファイルのファイル名
% timestep : numeric
%   タイムステップ
% u : matrix
%   可視化場
% maxrange : スカラー
%   コンターの最大値
% minrange : スカラー
%   コンターの最小値
% fignum : スカラー
% 　何番目の描写ウィンドウに書き込むか


% グローバル変数呼び出し
global dt

figure(fignum);
imagesc(u)
view(0, 90);%視点の設定
title(['time = ', num2str(timestep * dt, '%.3f')]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
axis equal; axis tight; axis on;
colorbar
caxis([maxrange minrange])
frame = getframe(fignum);
im = frame2im(frame);
[imind, cm] = rgb2ind(im, 256);
if timestep == 1
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'Loopcount', inf);
elseif rem(timestep, 10) == 0
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'WriteMode', 'append');
end

end

function[] = vis_vector(filename, timestep, u, v, fignum)

% グローバル変数呼び出し
global dt nx ny

figure(fignum);
quiver(flipud(rot90(u)),flipud(rot90(v)),'r')
title(['time = ', num2str(timestep * dt, '%.3f')]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
axis equal; axis tight; axis on;
xlim([0 nx]);
ylim([0 ny]);
frame = getframe(fignum);
im = frame2im(frame);
[imind, cm] = rgb2ind(im, 256);
if timestep == 1
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'Loopcount', inf);
elseif rem(timestep, 10) == 0
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'WriteMode', 'append');
end

end

function[uu, vv] = VelocityInterpolate(u, v, uu, vv)

% グローバル変数呼び出し
global nx ny

for i = 1 : nx + 2
    for j = 1 : ny + 2
        if i == 1 % 流入口ならば
            uu(i, j) = 0.5 * (3 * u(i, j) - u(i + 1, j));
        elseif  i == nx + 2 %流出口ならば
            uu(i, j) = 0.5 * (3 * u(i - 1, j) - u(i - 2, j));
        else% 内部領域
            uu(i, j) = 0.5 * (u(i, j) + u(i - 1, j));%
        end
    end
end
for i = 1 : nx + 2
    for j = 1 : ny + 2
        if j == 1
            vv(i, j) = 0.5 * (3 * v(i, j) - v(i, j + 1));%後退差分近似
        elseif j == ny + 2
            vv(i, j) = 0.5 * (3 * v(i, j - 1) - v(i, j - 2));%後退差分近似
        else
            vv(i, j) = 0.5 * (v(i, j) + v(i, j - 1));%前進差分近似
        end
    end
end

end

function[div, divup] = CheackContinuityEquation(up, vp, divup)

% グローバル変数呼び出し
global nx ny ddx ddy

ic = 0;
div = 0;
for j = 2 : ny + 1
    for i = 2 : nx + 1
        divup(i, j) = ddx * (up(i, j) - up(i - 1, j)) + ddy * (vp(i, j) - vp(i,j - 1));
        ic = ic + 1;
        div = div + divup(i, j)^2;
    end
end

end

function[u, v, p] = ModifyVP(up, vp, u, v, p, phi)

% グローバル変数呼び出し
global nx ny ddx ddy dt

for j = 2 : ny + 1
    for i = 2 : nx
        u(i, j) = up(i, j) - dt * ddx * (phi(i + 1, j)-phi(i, j));%式２７
    end
end
for j = 2 : ny
    for i = 2 : nx + 1
        v(i, j) = vp(i, j) - dt * ddy * (phi(i, j + 1) - phi(i, j));%式２９
    end
end
for j = 2 : ny + 1
    for i = 2 : nx + 1
        p(i, j) = p(i, j) + phi(i, j);%式３１
    end
end

end

function[phi] = PoissonSolver(alpha, phi, eps, maxitr, divup, nx, ny, ddt, ddx2, ddy2)

% グローバル変数使うと遅くなるから使わない。

for iter = 1 : maxitr% SOR法により圧力補正値を求める。
    error = 0;
    for j = 2 : ny + 1
        for i = 2 : nx + 1
            rhs = ddt * divup(i, j);%式２５右辺
            resid = ddx2 * (phi(i - 1,j) - 2 * phi(i, j) + phi(i + 1, j))...
                + ddy2 * (phi(i, j - 1) - 2 * phi(i, j)+phi(i, j + 1)) - rhs;
            dphi = alpha * resid / (2 * (ddx2 + ddy2));
            error = max(abs(dphi), error);
            phi(i, j) = phi(i, j) + dphi;%式２５をphi(i,j)についてまとめSOR法の形にしたもの
        end
    end
    
    % 境界条件の設定
    phi(1, 2 : ny + 1) = phi(2, 2 : ny + 1);%東側での圧力勾配０。
    phi(nx + 2, 2 : ny + 1) = 0;%西側境界条件
    phi(2 : nx + 1, 1) = phi(2 : nx + 1, 2);%南側境界条件
    phi(2 : nx + 1, ny + 2) = phi(2 : nx + 1, ny + 1); %北側境界条件
    
    if error < eps % 収束条件が満たされたらループを抜ける。
        break
    end
    
    if iter >= maxitr
        disp('最大反復回数に達しました。収束条件を満たしていません。');
    end
end

end

function[up] = ProvisionalVelocityU(u, v, p, up)

% グローバル変数呼び出し
global nx ny ddx ddy ddx2 ddy2 re dt

for j = 2 : ny + 1
    for i = 2 : nx % temporary u-velocity
        %u（ij）中心で計算
        %移流項の離散化
        cnvu = ddx * ((u(i + 1, j) + u(i, j))^2 - (u(i - 1, j) + u(i, j))^2) / 4 ...
            + ddy * ((u(i, j + 1) + u(i, j)) * (v(i + 1, j)+v(i, j))...
            -(u(i, j) + u(i, j - 1)) * (v(i, j - 1) + v(i + 1, j - 1))) / 4;
        fij = - ddx * (p(i + 1, j) - p(i, j)) - cnvu...
            + ddx2 * (u(i + 1, j) - 2 * u(i, j) + u(i - 1, j)) / re...
            + ddy2 * (u(i, j + 1) - 2 * u(i, j) + u(i, j - 1))/ re;
        up(i, j) = u(i, j) + dt * fij;
    end
end

end

function[vp] = ProvisionalVelocityV(u, v, p, vp)

% グローバル変数呼び出し
global nx ny ddx ddy ddx2 ddy2 re dt

for j = 2 : ny
    for i = 2 : nx + 1% temporary v-velocity
        % v（ij）中心でd計算
        % 移流項の離散化
        cnvv = ddx * ((u(i, j + 1) + u(i, j)) * (v(i + 1, j)+v(i, j))...
            - (u(i - 1, j + 1) + u(i - 1, j)) * (v(i - 1, j)+v(i, j))) / 4 ...
            + ddy * ((v(i, j + 1) + v(i, j))^2 - (v(i, j) + v(i, j - 1))^2) / 4;
        gij = - ddy * (p(i, j + 1) - p(i, j)) - cnvv...
            + ddx2 * (v(i + 1, j) - 2 * v(i, j) + v(i - 1, j)) / re...
            + ddy2 * (v(i,j + 1) - 2 * v(i, j) + v(i, j - 1)) / re;
        vp(i, j) = v(i, j) + dt * gij;
    end
end

end

function[u] = BoundaryConditionU(u, ue, uw, us, un)

% グローバル変数呼び出し
global nx ny

u(nx + 1, 1 : ny + 1) = u(nx,1 : ny + 1);% 速度勾配０
u(1, 1 : ny + 1) = uw; % 西側（流入側）境界条件
u(1 : nx + 1, 1) = us; % 南側境界条件
u(1 : nx + 1, ny + 2) = un; % 北側境界条件

end

function[v] = BoundaryConditionV(v, ve, vw, vs, vn)

% グローバル変数呼び出し
global nx ny

v(2 : nx + 1, 1) = vs;% 南側境界条件
v(2 : nx + 1, ny + 1) = vn;% 北側境界条件
v(1, 2 : ny) = vw;% 西側境界条件。端点は西側には含めず、東、南と考える。
v(nx + 2, 2 : ny) = v(nx + 1, 2 : ny);% 東側の速度勾配０

end

function[p] = BoundaryConditionP(p)

% グローバル変数呼び出し
global nx ny

p(nx + 2, 1 : ny + 1) = 0;% 東側（流出側）境界条件 圧力０
p(1, 1 : ny + 1) = p(2, 1 : ny + 1);% 西側（流入側）境界条件
p(1 : nx + 1, 1) = p(1 : nx + 1, 2);% 南側境界条件
p(1 : nx + 1, ny + 2) = p(1 : nx + 1, ny + 1);% 北側境界条件

end

function[object] = DefineObjectArea(object, center)

% グローバル変数呼び出し
global nx ny dx dy

% 障害物領域の定義
for i = 1 : nx + 2
    for j = 1 : ny + 2
        r = sqrt(((i - center(1)) * dx)^2 + ((j - center(2)) * dy)^2);%中心から格子点までの距離
        if r < 2.5 * dx
            object(i, j) = 1;% 障害物の位置を1とする。
        end
    end
end
% 障害物境界領域の抽出
[row1, col1] = find(object > 0);
for i = 1: size(row1)
    if object(row1(i) - 1, col1(i)) == 0 || object(row1(i), col1(i)-1)==0 ||...
            object(row1(i)+1, col1(i)) == 0 || object(row1(i),col1(i)+1) == 0
        object(row1(i), col1(i)) = 2;%角柱の境界を２とする。
    end
end

end
