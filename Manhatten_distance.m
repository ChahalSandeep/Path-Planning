%function prepared_world_map= constructing_worldmap_obstacle()
world_map=zeros(10,10);

%goal
world_map(8, 7)=3;

%obstacle one 
world_map(size(world_map,1)-2 : size(world_map,1), 3:4)=1;   % obstacle
world_map(size(world_map,1) - 3, 3:4) = 2;                  % extended obstacle  
world_map(size(world_map,1) - 3:size(world_map,1), 2) = 2;  % extended obstacle
world_map(size(world_map,1) - 3:size(world_map,1), 5) = 2;  % extended obstacle

%Obstacle 2
world_map(2:3, 7:8) = 1; % obstacle
world_map(1, 7:8) = 2;  % extended obstacle
world_map(4, 7:8) = 2;  % extended obstacle
world_map(1:4, 6) = 2;  % extended obstacle
world_map(1:4, 9) = 2;  % extended obstacle

occuGrid=world_map;
%fucntion 
ix=1;    %initial x location in the grid
iy=10;    %initial x location in the grid  
heading=0;
%k=find(world_map==3);
[row,col] = find(world_map==3);
gx=col;
gy=row;
%function  manhatten_grid = manhatten(gx,gy,occuGrid)
manhatten_grid = zeros(size(occuGrid, 1), size(occuGrid, 2));
manhatten_grid(occuGrid==1)=Inf; % setting manhatten grif to infinite at places og obstacle and ectened obstacle i.e. 1 and 2
manhatten_grid(occuGrid==2)=Inf;
parse_next=0;
%logical(sum(manhatten_grid == 0)~= 1)
a=numel(manhatten_grid);  
while (a-nnz(manhatten_grid))~=1
    if parse_next ==0
        
        % function manhatten_grid = mark_cells(y,x,parse_next,manhatten_grid,occuGrid);
        distance = parse_next;
        if (row ~= size(manhatten_grid, 1))
            if (manhatten_grid(row + 1, col) == 0) && (occuGrid(row + 1, col) == 0)
                manhatten_grid(row + 1, col) = distance + 1;
            end
        end
        
        if (row ~= 0)
            if ((manhatten_grid(row - 1, col) == 0) && (occuGrid(row - 1, col) == 0))
                manhatten_grid(row - 1, col) = distance + 1;
            end
        end
        if (col ~= size(manhatten_grid, 2))
            if ((manhatten_grid(row, col + 1) == 0) && (occuGrid(row, col + 1) == 0))
                manhatten_grid(row, col + 1) = distance + 1;
            end
        end
        if (col ~= 0)
            if ((manhatten_grid(row, col - 1) == 0) && (occuGrid(row, col - 1) == 0))
                manhatten_grid(row, col - 1) = distance + 1;
            end
        end
        %return manhatten_grid
        
        parse_next = parse_next+1;
    end
    [row,col]=find(manhatten_grid==parse_next);
    next_indices=[row,col];
    for index = 1:length(next_indices)
        %fucntion manhatten_grid = mark_cells(y,x,parse_next,manhatten_grid,occuGrid);
        row=next_indices(index,1);
        col=next_indices(index,2);
        distance = parse_next;
        if (row ~= size(manhatten_grid, 1))
            if (manhatten_grid(row + 1, col) == 0) && (occuGrid(row + 1, col) == 0)
                manhatten_grid(row + 1, col) = distance + 1;
            end
        end
        
        if (row ~= 1)
            if ((manhatten_grid(row - 1, col) == 0) && (occuGrid(row - 1, col) == 0))
                manhatten_grid(row - 1, col) = distance + 1;
            end
        end
        if (col ~= size(manhatten_grid, 2))
            if ((manhatten_grid(row, col + 1) == 0) && (occuGrid(row, col + 1) == 0))
                manhatten_grid(row, col + 1) = distance + 1;
            end
        end
        if (col ~= 1)
            if ((manhatten_grid(row, col - 1) == 0) && (occuGrid(row, col - 1) == 0))
                manhatten_grid(row, col - 1) = distance + 1;
            end
        end
        %return manhatten_grid
    end
    parse_next = parse_next+1;
    
end
disp ('Occupancy Grid')
disp(occuGrid)
disp('Manhatten Grid')
disp(manhatten_grid)
path = [];
path{end+1} = [iy,ix];

%lowest selection function
row=iy;
col=ix;
minimum=Inf;
min_index=[1,1];
if row~=size(manhatten_grid, 1)
    minimum=manhatten_grid(row+1,col);
    min_index=[row+1,col];
end
if row~=1
    if minimum > manhatten_grid(row-1,col)
        minimum = manhatten_grid(row-1,col);
        min_index=[row-1,col];
    end
end
if col~=size(manhatten_grid, 2)
    if minimum>manhatten_grid(row,col+1)
        minimum = manhatten_grid(row,col+1);
        min_index=[row,col+1];
    end
end
if col~=1
    if minimum>manhatten_grid(row,col-1)
        minimum = manhatten_grid(row,col-1);
        min_index=[row,col-1];
    end
end
%return minimum and min_index
distance=minimum;
path{end+1} = min_index;
index=min_index;
while distance~=0
    %function
    row=index(1);
    col=index(2);
    minimum=Inf;
    min_index=[1,1];
    if row~=size(manhatten_grid, 1)
        minimum=manhatten_grid(row+1,col);
        min_index=[row+1,col];
    end
    if row~=1
        if minimum > manhatten_grid(row-1,col)
            minimum = manhatten_grid(row-1,col);
            min_index=[row-1,col];
        end
    end
    if col~=size(manhatten_grid, 2)
        if minimum>manhatten_grid(row,col+1)
            minimum = manhatten_grid(row,col+1);
            min_index=[row,col+1];
        end
    end
    if col~=1
        if minimum>manhatten_grid(row,col-1)
            minimum = manhatten_grid(row,col-1);
            min_index=[row,col-1];
        end
    end
   distance=minimum;
   index=min_index;
   path{end+1} = [min_index];
end
disp('Way points')
%disp(path)
path_matrix=zeros(size(manhatten_grid));
for node=1:length(path)
    temp=path(node);
    a=temp{1}(1);
    b=temp{1}(2);
    disp([a,b])
    path_matrix(a,b)=1;
end
disp('PATH')
disp(path_matrix)