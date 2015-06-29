#A julia file to organize coarsegraining of the DFAX series
#rename all residues of .pdb file
function parseDict(filename)
r = open(filename);
resDict = (String => String)[];
ls = readlines(r);
for line in ls

	l1 = split(line);
	if length(l1) > 0
	#println(l1)
	key = l1[1];
	val = l1[2];
	#if val[end]=="r"
	#val = join([" ",val]);
	#end
	resDict[key] = val;
	end
end
#println(resDict)
return resDict
end

function renameRes()
resDict = parseDict("resDict.txt");
cadict = parseDict("CAdict.txt");
println(resDict);
f = open("testDFAI.pdb");
o = open("testDFAImod.pdb","w");
for line in eachline(f)
	#println(line);
	if contains(line,"ATOM")
		splitline = split(line);
		if haskey(resDict,splitline[4])
			key = splitline[4];
			val = resDict[splitline[4]];
		else
			val = splitline[4];
		end

		if val[end] == 'r'
			key = join([key," "]);
		end
		if length(key) < length(val)
			key = join([repeat(" ",length(val)-length(key)),key]);
		end
		outline = replace(line,key,val);
		#println(val);
		#println(splitline[3]);
		if ((val == "COP") && ((splitline[3] == "C") || (splitline[3]) == "O"))			
			#println("dealing with COA");
			#take care of the special COA residue
			val = "COA";
			outline = replace(outline,"COP","COA");
			if splitline[3] == "C"
				outline = replace(outline,"C ","CA");
			end
		end
		if haskey(cadict,val)
		#take care of incorrect assignment of CA atoms
		#println(val)
		#println(splitline[3])
		#println(cadict[val])
			if splitline[3] == cadict[val]
				#println("dealing with CA atoms");
				outline = replace(outline,cadict[val],"CA");
				
			end
		end
	#	println(outline);
		write(o,outline);
	else
		outline = line;
		#println(outline);
		write(o,line);
	end
end
#println("Go edit formatting on testDFAImod.pdb.  Still working on this stupid bug.");
end
renameRes();
#lines = readlines(f);


#run martinize to get the SC correct
#run(`bash testmart.sh`);
#run(`python martinize_hacked.py -f testDFAImod.pdb -x CG_DFAI.pdb -n CG.ndx -nmap mapping.ndx -nt -ss ss.dat -o CG_DFAI.top -ff martini22`);

println("Currently, need to run martinize_hacked.py by hand to get this thing to work.");



#fix bonds/angles/dihedrals in the .itp file
