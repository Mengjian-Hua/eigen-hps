function out = isInitialized(S)
%ISINITIALIZED   Check to see if a SURFACEOP has been initialized.

out = ~isempty(S.patches) && (~isempty(S.patches{1}.Iu_part) );

end
