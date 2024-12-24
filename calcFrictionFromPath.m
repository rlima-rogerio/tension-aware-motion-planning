function [totalFriction, tensionPerSegment] = calcFrictionFromPath(ge, path_, sizeVisTable, hashtable)

totalFriction = 0;
tensionPerSegment = [];

% Includes the cost of the second edge of the virtual start (0,1) -> (1,2)
pt1 = 0;
pt2 = path_(1);
pt3 = path_(2);
nodeIdx1 = ge.lookup(sprintf('%d', matrix2digit(pt1, pt2, sizeVisTable)));
nodeIdx2 = ge.lookup(sprintf('%d', matrix2digit(pt2, pt3, sizeVisTable)));
key = sprintf("%d.%d", nodeIdx1, nodeIdx2);
edgeID = hashtable(key);
cost = ge.cost(edgeID);
totalFriction = totalFriction + cost;
tensionPerSegment = [tensionPerSegment, totalFriction];

for ptItr = 2 : length(path_) - 1

    pt1 = path_(ptItr - 1);
    pt2 = path_(ptItr);
    pt3 = path_(ptItr + 1);

    nodeIdx1 = ge.lookup(sprintf('%d', matrix2digit(pt1, pt2, sizeVisTable)));
    nodeIdx2 = ge.lookup(sprintf('%d', matrix2digit(pt2, pt3, sizeVisTable)));

    key = sprintf("%d.%d", nodeIdx1, nodeIdx2);
    edgeID = hashtable(key);
    cost = ge.cost(edgeID);

    totalFriction = totalFriction + cost;
    tensionPerSegment = [tensionPerSegment, totalFriction];
end

% Includes the cost of the first edge of the virtual goal (n-1,n) -> (n,n+1)
pt1 = path_(end-1);
pt2 = path_(end);
pt3 = sizeVisTable + 1;
nodeIdx1 = ge.lookup(sprintf('%d', matrix2digit(pt1, pt2, sizeVisTable)));
nodeIdx2 = ge.lookup(sprintf('%d', matrix2digit(pt2, pt3, sizeVisTable)));
key = sprintf("%d.%d", nodeIdx1, nodeIdx2);
edgeID = hashtable(key);
cost = ge.cost(edgeID);
totalFriction = totalFriction + cost;
tensionPerSegment = [tensionPerSegment, totalFriction];