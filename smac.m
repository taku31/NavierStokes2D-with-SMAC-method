% Solving 2D navier stokes equation with SMAC method
% Copyright (C) 2019  T.Nakabayashi
% Released under the MIT license http://opensource.org/licenses/mit-license.php

% 初期化
clear all;

% グローバル変数宣言
global dt ddt nx ny dx dy ddx ddx2 ddy ddy2 re

% パラメーター
n = 60;% 格子数
nx = 2 * n;% x方向格子数
ny = 1 * n;% % ｙ方向格子数
loop = 3000;% ステップ数
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
omega = zeros(nx + 2, ny + 2);% 渦度
uu = zeros(nx+2,ny+2);
vv = zeros(nx+2,ny+2);
pp = zeros(nx+2,ny+2);
omega2 = zeros(nx + 2, ny + 2);

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
u(:,ny+2) = 1;% 全領域を1にすることで、初タイムステップでできる限り連続の方程式を満たすようにしている？

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
center = round([(nx + 2) / 5, (ny + 2) / 2]);% 中心節点番号。整数指定しないと、非対称になる。
r = 6;% 半径（節点個数単位）
object = DefineObjectArea(object, center, r);

% 粒子位置の生成
%rng(1);% 乱数シードの固定
p_num = 30;% 1タイムステップで発生させる粒子数。
dt_p = 50 * dt;% 粒子発生間隔
tmax_p = 200;% 粒子発生完了時間
p_num_total = p_num * (fix(tmax_p / dt_p) + 1);% 合計発生粒子数。初期配置粒子分1を足しておく。
px = NaN(p_num_total, 1);% 粒子x座標。NaN要素として確保すると、描写時に粒子がない要素を無視できる。
py = NaN(p_num_total, 1);% 粒子y座標。

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
    eps = 10^(- 10);
    maxitr = nx * ny * 100;% 反復回数。収束させるためにはこのぐらい必要。
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
    
    % 渦度の計算
    omega = CalVorticity(u, v, omega);
    
    % 圧力格子位置での速度を求める。
    [uu, vv, pp, omega2] = interpolation(u, v, p, omega, uu, vv, pp, omega2, object);
    
    % 粒子の生成
    [px, py] = GenerateParticles(px, py, p_num, tmax_p, dt_p, ita);
    
    % 粒子位置の更新
    [px, py] = FlowParticles(px, py, uu, vv);
    
    % 結果の描画
    %vis_contour('u.gif', ita, uu, 0, 1.5, 1)
    %vis_contour('v.gif', ita, vv, -0.5, 0.5, 2)
    %vis_contour('p.gif', ita, pp, -0.5, 0.5, 3)
    %vis_contour('vorticity.gif', ita, omega2, -2.5, 2.5, 4)
    %vis_contour('div.gif', ita, divup, -0.0001, 0.0001, 5)
    %vis_vector('vec.gif', ita, uu, vv, 6)
    vis_particles('particle.gif', ita, px, py, object, 7)
    
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
global dt nx ny

figure(fignum);
u = flipud(rot90(u));
h = imagesc(u);
view(0, 270);%視点の設定
h.AlphaData = isfinite(u); % NaNやInfを透明にする
title(['time = ', num2str(timestep * dt, '%.3f')]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
axis equal; axis tight; axis on;
colorbar('southoutside')

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

function[] = vis_particles(filename, timestep, px, py, object, fignum)

% グローバル変数呼び出し
global dt nx ny

figure(fignum);
clf(fignum,'reset')
hold on;
object2 = flipud(rot90(object));
object2(object2 == 2) = 0;
imagesc(object2);
c = gray;
c = flipud(c);
colormap(c);
scatter(px, py)
title(['time = ', num2str(timestep * dt, '%.3f')]);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 16);
axis equal; axis tight; axis on;
xlim([0 nx]);
ylim([0 ny]);
hold off;

frame = getframe(fignum);
im = frame2im(frame);
[imind, cm] = rgb2ind(im, 256);
if timestep == 1
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'Loopcount', inf);
elseif rem(timestep, 10) == 0
    imwrite(imind, cm, filename, 'gif', 'DelayTime', 0.001, 'WriteMode', 'append');
end

end

function[uu, vv, pp, omega2] = interpolation(u, v, p, omega, uu, vv, pp, omega2, object)

% グローバル変数呼び出し
global nx ny


for i = 1 : nx + 2
    for j = 1 : ny + 2
        if i == 1
            uu(i, j) = 0.5 * (3 * u(i, j) - u(i + 1, j));
        elseif  i == nx + 2
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

for i = 1 : nx + 2
    for j = 1 : ny + 2
        if i == 1
            omega2(i, j) = 0.5 * (3 * omega(i, j) - omega(i + 1, j));
        elseif j == 1
            omega2(i, j) = 0.5 * (3 * omega(i, j) - omega(i, j + 1));%後退差分近似
        elseif  i == nx + 2
            omega2(i, j) = 0.5 * (3 * omega(i - 1, j) - omega(i - 2, j));
        elseif j == ny + 2
            omega2(i, j) = 0.5 * (3 * omega(i, j - 1) - omega(i, j - 2));%後退差分近似
        else
            omega2(i, j) = 0.5 * (omega(i, j) + omega(i, j - 1));%前進差分近似
        end
    end
end

pp = p;%圧力は補間の必要なし。障害物処理のみする。

%障害物領域はNANにする。
for i = 1 : nx + 2
    for j = 1 : ny + 2
        if object(i, j) == 1
            uu(i, j) = NaN;
            vv(i, j) = NaN;
            pp(i, j) = NaN;
            omega2(i, j) = NaN;
        end
    end
end

end

function[div, divup] = CheackContinuityEquation(up, vp, divup)

% グローバル変数呼び出し
global nx ny ddx ddy

ic = 0;
divsum = 0;

for j = 2 : ny + 1
    for i = 2 : nx + 1
        divup(i, j) = ddx * (up(i, j) - up(i - 1, j)) + ddy * (vp(i, j) - vp(i,j - 1));
        ic = ic + 1;
        divsum = divsum + divup(i, j)^2;
    end
end

div = sqrt(divsum/ic);

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

function[object] = DefineObjectArea(object, center, r)

% グローバル変数呼び出し
global nx ny dx dy

% 障害物領域の定義(円形)
for i = 2 : nx + 1
    for j = 2 : ny + 1
        radius = sqrt(((i - center(1))* dx)^2 + ((j - center(2)) * dy)^2 );%中心から格子点までの距離
        if radius < r * dx
            object(i, j) = 1;% 障害物の位置を1とする。
        end
    end
end

% % 障害物領域の定義(四角)
% for i = 1 : nx + 2
%     for j = 1 : ny + 2
%         if (i - 2) * dx > center(1) - r && (i - 2) * dx < center(1) + r && (j - 2) * dy > center(2) - r && (j - 2) * dy < center(2) + r
%             object(i, j) = 1;% 障害物の位置を1とする。
%         end
%     end
% end

% 障害物境界領域の抽出
[row1, col1] = find(object > 0);
for i = 1: size(row1)
    if object(row1(i) - 1, col1(i)) == 0 || object(row1(i), col1(i)-1)==0 ||...
            object(row1(i)+1, col1(i)) == 0 || object(row1(i),col1(i)+1) == 0
        object(row1(i), col1(i)) = 2;%障害物の境界を２とする。
    end
end

end

function[omega] = CalVorticity(u, v, omega)

% グローバル変数呼び出し
global nx ny ddx ddy

for j = 2 : ny + 1
    for i = 2 : nx + 1
        omega(i, j) = ddx * (v(i, j) - v(i - 1, j)) - ddy * (u(i, j) - u(i,j - 1));
    end
end

end


function[px, py] = FlowParticles(px, py, uu, vv)

% グローバル変数呼び出し
global  dt nx ny

for i = 1 : size(px, 1)
    
    if ~isnan(px(i))% 粒子が未存在の場合は処理しない。
        
        px_int = fix(px(i));% 粒子x座標の0に近い側の整数
        py_int = fix(py(i));% 粒子y座標の0に近い側の整数
        
        px_s= px(i) - px_int; % 粒子x座標の小数部分
        py_s= py(i) - py_int; % 粒子y座標の小数部分
        
        % 簡易的な補間計算より移動量を求める。
        spdx = (px_s * uu(px_int, py_int) + px_s * uu(px_int, py_int + 1)...
            + (1 - px_s) * uu(px_int + 1, py_int ) + (1 - px_s) * uu(px_int + 1, py_int + 1)) * dt;
        spdy = (py_s * vv(px_int, py_int) + py_s * vv(px_int, py_int + 1)...
            + (1 - py_s) * vv(px_int + 1, py_int ) + (1 - py_s) * vv(px_int + 1, py_int + 1)) * dt;
        
        % もしNANならば粒子は動かさない。
        if isnan(spdx)
            spdx = 0;
        end
        if isnan(spdy)% elseif使うとうまく判定されなかったのでifで。
            spdy = 0;
        end
        
        % 粒子位置の更新
        px(i) = px(i) + spdx;
        py(i) = py(i) + spdy;
        
        % 領域外に出た粒子は削除する。
        if px(i) >= nx + 2 || px(i) < 1 || py(i) >= ny + 2 || py(i) < 1
            px(i) = NaN;% NaNにしておくとscatter描写時に無視される。
            py(i) = NaN;
        end
        
    end
    
end

end

function[px, py] = GenerateParticles(px, py, p_num, tmax_p, dt_p, ita)

% グローバル変数呼び出し
global  nx ny dt

if ita * dt <= tmax_p
    if  rem(ita * dt, dt_p) == 0 || ita == 1
        
        pc = fix(1 + (ita * dt) / dt_p);% 粒子の発生回数
        
        %粒子のx座標を生成。
        px((pc - 1) * p_num + 1 : (pc - 1) * p_num + p_num) = 1;
        
        %粒子のy座標を生成。
        py_gene = 1 : (ny + 2) / (p_num + 1) : ny + 2;
        py((pc - 1) * p_num + 1 : (pc - 1) * p_num + p_num) = py_gene(1 : p_num);
        
    end
end

end
