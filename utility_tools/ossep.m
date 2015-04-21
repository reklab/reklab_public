function sep=ossep(computer);

switch computer,
    case 'PCWIN'
        sep= '\';
    case {'SOL2','GLNX86'}
        sep= '/';
    otherwise
        error (['Computer type not supported:' computer] );
end