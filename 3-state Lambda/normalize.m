function v=normalize(x)
    md=norm(x);
    v=simplify(x./md,'Steps',100);
end