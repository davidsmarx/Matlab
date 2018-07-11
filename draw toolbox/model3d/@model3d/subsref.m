function b = subsref(S, s)

switch s(1).type,
    case '()'
        b = S(s(1).subs{1});
        
    case '{}'
        
    case '.'
        
        stmp = struct(S);
        b = stmp.(s(1).subs);
        
end % switch s.type

if length(s) > 1,
    s2 = s(2:end);
    b = subsref(b,s2);
end

end % function
