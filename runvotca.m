function v = runvotca(votcainputname,mdtop,mdtrj,outxml)
    comm = ['(cat ' votcainputname ') | csg_boltzmann --top ' mdtop ' --trj ' mdtrj ' --cg "' outxml ';water.xml"'];
    system(comm);