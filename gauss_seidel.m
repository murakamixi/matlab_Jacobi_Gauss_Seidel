function [outputX,outputRes] = gauss_seidel(A, b, torelance)
% gauss_seidel法の計算
% A:n次正則な行列
% b:定数
% torelance:許容誤差

% 初期値（適当な解）x0と残差を設定
x = [0 0 0];
residual = 1e20;

D = transpose(diag(A)); % D := 対角行列
R = A - diag(D); % R := L + U

iteration = 0;

disp('----------Start iteration----------');
res = [];

% 反復計算の結果、誤差が許容誤差（torelance）よりも小さくなったら終了
while residual > torelance

    % 収束判定用に1ステップ前の解を保持
    x0 = x;
    for i=1:length(x)
        % R(i, :)でi番目の行ごと取り出す
        x(i) = (b(i) - transpose(R(i, :)*transpose(x))) ./ D(i);
    end

    % xがbでx0がAxのつもりで、どれだけ解が欲しいターゲットの値と乖離しているかを計算
    denominator = sqrt((x - x0).^2);
    numerator = sqrt(x.^2);
    residual = sum(denominator) / sum(numerator);
    res = [res residual];

    % 100回ごとに誤差を出力
    if rem(iteration, 100) == 0
        fprintf('Iteration = %d \n', iteration);
        fprintf('Residual = %d \n', residual);
        fprintf('x = %d \n', x);
    end

    % 更新した結果を変数に格納して次の処理に持ち越す。
    iteration = iteration + 1;
end

% 終了時の結果を出力
disp('----------End iteration----------')
fprintf('Iteration = %d \n', iteration);
fprintf('Residual = %d \n', residual);
fprintf('x = %d \n', x0);

outputX = x0;
outputRes = res;
end