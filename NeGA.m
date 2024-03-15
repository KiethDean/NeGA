close;
clear;
clc;

%% 可调参数
% 元胞自动机模拟生态系统的重要参数
n = 100;                      % 元胞（方）矩阵大小（100即为10000个格子）
UL = [n 1:n-1];               % 设置上/左边界
DR = [2:n 1];                 % 设置下/右边界
%veg = randi([0, 5], n, n);   % 初始化为随机矩阵
veg = zeros(n);               % 初始化为全零矩阵
Disaster = 0.00;              % 灾难发生频率（越接近1，灾害发生频率越大）

Availability_factor = 0.30;   % 资源可用性（随季节波动）对性别比例的影响（0到0.5，数值越大影响越强）
                              % （ps：这里后面可能还需要做工作，把资源可用性与Availability_factor关联起来，映射一个0到0.5的区间。）

GrowFactor_impact = 15;       % 性别比例对生长因子的影响（默认值为10）

     %（ps：各物种受周围邻居的影响的量化在后面，后面的部分谨慎改参）

% 物种的数量和标识符、抗灾能力
species_anti_disaster = [2, 1.1, 1.0, 1.05, 1];      % 5个物种的抗灾能力（数值越大，抗灾能力越强，0代表没有抗灾能力）
                %  生产者↑           ↑最高级消费者
                %  下方顺序同理
speciesCounts = [3000, 1000, 700, 1000, 1200]; % 5个物种，每个物种的个体数量
%%
speciesIDs = 1:5; % 物种的标识符

% 随机分配每个物种的位置
for i = 1:length(speciesIDs)
    speciesID = speciesIDs(i);
    count = speciesCounts(i);
    
    % 对于当前物种，随机选择位置并赋值
    for j = 1:count
        % 随机选择一个位置
        while true
            row = randi(n);
            col = randi(n);
            % 如果这个位置已经被占用，就重新选择
            if veg(row, col) == 0
                veg(row, col) = speciesID;
                break; % 跳出循环，继续下一个个体
            end
        end
    end
end

% 此时matrix矩阵中包含了随机分布的6个物种
% 初始化图像显示
f1 = figure(1);
movegui(f1,"northwest")

imh = imagesc(veg);
color = jet(6);
colormap(color); % 使用颜色图以区分6种状态
colorbar;
axis equal tight;


f2 = figure(2); clf;
movegui(f2,"northeast")
hold on; % 保持图表，以便添加新的绘图线

f3 = figure(3); clf;
hold on; % 保持图表，以便添加新的绘图线

m=annotation('textbox',[0.1,0.1,0.1,0.1],'LineStyle','-','LineWidth',1,'String','123');

for i = 1:300
    %性别比例
    %Ratio(i) = (1./(1+exp(1-0.8 + 0.2 * sin(i * pi))));

    %考虑性别比例后的生长因子
    GrowFactor(i) = 1 - GrowFactor_impact * ((1./(1+exp(1-(1-Availability_factor) + Availability_factor * sin(i * pi * 0.2)))) - 0.5)^2; 
    %不考虑性别比例后的生长因子
    %GrowFactor(i) = 1;
    
    % 更新规则示例（需要根据实际生态系统逻辑调整）
    new_veg = veg;
    for row = 1:n
        for col = 1:n
            neighbor = getRandomNeighbor(veg, row, col);
            switch veg(row, col)
                case 0 % 空地
                    if rand < 0.8 % 生产者自然生长概率
                        new_veg(row, col) = 1;
                    end

                case 1 % 生产者
                    switch neighbor
                        case 2 || 3 || 4
                            random_number = rand;
                            if random_number < 0.5
                                new_veg(row, col) = 0;
                            elseif random_number < 0.7 * GrowFactor(i)
                                new_veg(row, col) = 2;
                            elseif random_number < 0.85 * GrowFactor(i)
                                new_veg(row, col) = 3;
                            elseif random_number < 1.0
                                new_veg(row, col) = 4;
                            end

                        case 5
                            if rand < 0.2
                                new_veg(row, col) = 5;
                            end

                        otherwise
                            random_number = rand;
                            if random_number < 0.15
                                new_veg(row, col) = 0;
                            elseif random_number < 0.26
                                new_veg(row, col) = 2;
                            end
                    end
                    if new_veg(row, col) ~= 0
                        die_number = species_anti_disaster(new_veg(row, col)) * rand;
                        if die_number < Disaster
                            new_veg(row, col) = 0;
                        end
                    end


                case 2 % 七鳃鳗
                    switch neighbor
                        case 3
                            if rand < 0.5
                                new_veg(row, col) = 3;
                            end
                        case 4
                            random_number = rand;
                            if random_number < 0.1
                                new_veg(row, col) = 0;
                            elseif random_number < 0.2
                                new_veg(row, col) = 1;
                            elseif random_number < 0.5
                                new_veg(row, col) = 2;
                            elseif random_number < 0.6
                                new_veg(row, col) = 3;
                            elseif random_number < 0.95
                                new_veg(row, col) = 4;    
                            elseif random_number < 1.0
                                new_veg(row, col) = 5;
                            end
                        case 5
                            random_number = rand;
                            if random_number < 0.4
                                new_veg(row, col) = 0;
                            elseif random_number < 0.55
                                new_veg(row, col) = 3;
                            elseif random_number < 0.58
                                new_veg(row, col) = 4;
                            elseif random_number < 1.0
                                new_veg(row, col) = 5;
                            end
                        
                        otherwise
                            random_number = rand;
                            if random_number < 0.15
                                new_veg(row, col) = 0;
                            elseif random_number < 0.3
                                new_veg(row, col) = 4;
                            elseif random_number < 0.315
                                new_veg(row, col) = 5;
                            end
                    end
                    if new_veg(row, col) ~= 0
                        die_number = species_anti_disaster(new_veg(row, col)) * rand;
                        if die_number < Disaster
                            new_veg(row, col) = 0;
                        end
                    end

                case 3 % 共生者
                    switch neighbor
                        case 2
                            if rand < 0.6
                                new_veg(row, col) = 2;
                            end
                        case 4
                            random_number = rand;
                            if random_number < 0.5
                                new_veg(row, col) = 0;
                            elseif random_number < 1.0
                                new_veg(row, col) = 3;
                            end
                        case 5
                            random_number = rand;
                            if random_number < 0.4
                                new_veg(row, col) = 0;
                            elseif random_number < 1.0
                                new_veg(row, col) = 5;
                            end

                        otherwise
                            random_number = rand;
                            if random_number < 0.2
                                new_veg(row, col) = 0;
                            elseif random_number < 0.21
                                new_veg(row, col) = 5;
                            end
                    end
                    if new_veg(row, col) ~= 0
                        die_number = species_anti_disaster(new_veg(row, col)) * rand;
                        if die_number < Disaster
                            new_veg(row, col) = 0;
                        end
                    end
                              
                case 4 % 竞争者
                    switch neighbor
                        case 2
                            random_number = rand;
                            if random_number < 0.1
                                new_veg(row, col) = 0;
                            elseif random_number < 0.2
                                new_veg(row, col) = 1;
                            elseif random_number < 0.9
                                new_veg(row, col) = 2;
                            elseif random_number < 1.0
                                new_veg(row, col) = 4;
                            end
                        case 3
                            random_number = rand;
                            if random_number < 0.3
                                new_veg(row, col) = 2;
                            elseif random_number < 0.7
                                new_veg(row, col) = 3;
                            elseif random_number < 1.0
                                new_veg(row, col) = 4;
                            end
                        case 4
                            if rand < 0.4
                                new_veg(row, col) = 0;
                            end
                        case 5
                            random_number = rand;
                            if random_number < 0.3
                                new_veg(row, col) = 0;
                            elseif random_number < 1.0
                                new_veg(row, col) = 5;
                            end

                        otherwise
                            random_number = rand;
                            if random_number < 0.1
                                new_veg(row, col) = 0;
                            elseif random_number < 0.101
                                new_veg(row, col) = 2;
                            elseif random_number < 0.1025
                                new_veg(row, col) = 5;
                            end
                    end
                    if new_veg(row, col) ~= 0
                        die_number = species_anti_disaster(new_veg(row, col)) * rand;
                        if die_number < Disaster
                            new_veg(row, col) = 0;
                        end
                    end

                case 5 % 最高级消费者
                    switch neighbor
                        case 2
                            if rand < 0.5
                                new_veg(row, col) = 2;
                            end
                        case 3
                            random_number = rand;
                            if random_number < 0.53
                                new_veg(row, col) = 3;
                            end
                        case 4
                            if rand < 0.5
                                new_veg(row, col) = 4;
                            end

                        otherwise
                            random_number = rand;
                            if random_number < 0.12
                                new_veg(row, col) = 0;
                            elseif random_number < 0.123
                                new_veg(row, col) = 2;
                            end
                    end
                    if new_veg(row, col) ~= 0
                        die_number = species_anti_disaster(new_veg(row, col)) * rand;
                        if die_number < Disaster
                            new_veg(row, col) = 0;
                        end
                    end
            end
        end
        
    end
    %刷新一次 分布
    veg = new_veg;
    %计算一次 生态系统里各个物种的数目 并汇总到eco中
    num_empty(i,:) = length(find(veg==0));
    eco(i,1) = num_empty(i,:);
    num_producer(i,:) = length(find(veg==1));
    eco(i,2) = num_producer(i,:);
    num_lamprey(i,:) = length(find(veg==2));
    eco(i,3) = num_lamprey(i,:);
    num_symbiont(i,:) = length(find(veg==3));
    eco(i,4) = num_symbiont(i,:);
    num_competitor(i,:) = length(find(veg==4));
    eco(i,5) = num_competitor(i,:);
    num_toppredator(i,:) = length(find(veg==5));
    eco(i,6) = num_toppredator(i,:);
    
    % 评价部分
    % 计算生物多样性指数(香农-威纳指数)
    H(i,:) = 0;
    P(i,1) = eco(i,2) / n^2; 
    P(i,2) = eco(i,3) / n^2; 
    P(i,3) = eco(i,4) / n^2; 
    P(i,4) = eco(i,5) / n^2;
    P(i,5) = eco(i,6) / n^2; 
    for j = 1:5
        H(i,:) = H(i,:) - P(i,j) * log(P(i,j));
    end
    figure(3);
    plot(H,'Color',color(1,:));
    drawnow;

    % Figure 1 生长因子绘图
    str1 = "生长因子";
    str2 = GrowFactor(i);
    str = [str1 str2];

    % 更新图像
    set(imh, 'CData', veg);
    figure(2);
    %delete(m);

    plot(num_empty,'Color',color(1,:),'Linewidth', 2);
 
    plot(num_producer,'Color',color(2,:),'Linewidth', 2);

    plot(num_lamprey,'Color',color(3,:),'Linewidth', 2);

    plot(num_symbiont,'Color',color(4,:),'Linewidth', 2);

    plot(num_competitor,'Color',color(5,:),'Linewidth', 2);

    plot(num_toppredator,'Color',color(6,:),'Linewidth', 2);
    f2l = legend(['Empty        ',num2str(num_empty(i))],['Producer     ',num2str(num_producer(i))],['Lamprey      ',num2str(num_lamprey(i))],['Symbiont     ',num2str(num_symbiont(i))],['Competitor   ',num2str(num_competitor(i))],['Top-Predator ',num2str(num_toppredator(i))]);
    set(f2l,'FontName','FixedWidth')
    title(['Time = ',num2str(0.1 * i),' Year(s)']);
    %m=annotation('textbox',[0.15,0.8,0.1,0.1],'LineStyle','-','LineWidth',1,'String',str);
    drawnow;
   
    % 暂停以便观察
    pause(0.01);
end

function neighbor = getRandomNeighbor(matrix, row, col)
    [rows, cols] = size(matrix);
    
    % 确定邻居的可能的行和列偏移
    rowOffsets = [-1, 0, 1];
    colOffsets = [-1, 0, 1];
    
    % 存储可能的邻居位置
    neighbors = [];
    
    % 遍历所有可能的偏移，找到有效的邻居位置
    for i = 1:length(rowOffsets)
        for j = 1:length(colOffsets)
            % 排除自身位置的偏移
            if rowOffsets(i) == 0 && colOffsets(j) == 0
                continue;
            end
            
            newRow = row + rowOffsets(i);
            newCol = col + colOffsets(j);
            
            % 检查新位置是否在矩阵边界内
            if newRow >= 1 && newRow <= rows && newCol >= 1 && newCol <= cols
                % 将有效的邻居位置加入列表
                neighbors = [neighbors; newRow, newCol];
            end
        end
    end
    
    % 如果没有有效的邻居，返回空
    if isempty(neighbors)
        neighbor = [];
        return;
    end
    
    % 从邻居列表中随机选取一个邻居
    numNeighbors = size(neighbors, 1);
    chosenIndex = randi(numNeighbors);
    chosenNeighborRow = neighbors(chosenIndex, 1);
    chosenNeighborCol = neighbors(chosenIndex, 2);
    
    % 返回选中的邻居的值
    neighbor = matrix(chosenNeighborRow, chosenNeighborCol);
end