function nbid = neighborhood(sid,ncol,p)

id1 = sid-ncol;
if id1<=0
    id1 = [];
end

id2 = sid+ncol;
if id2>p
    id2 = [];
end

if rem(sid,ncol)~=1
    id3 = sid-1;
else
    id3 = [];
end

if rem(sid,ncol)~=0
    id4 = sid+1;
else
    id4 = [];
end

nbid = [id1,id2,id3,id4];