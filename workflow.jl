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

resDict = parseDict("resDict.txt");
#-
println(resDict);
f = open("testDFAI.pdb");
o = open("testDFAImod.pdb","w");
for line in eachline(f)
	if contains(line,"ATOM")
		splitline = split(line);
		if haskey(resDict,splitline[4])
			key = splitline[4];
			val = resDict[splitline[4]];

		#if length(key) > length(val)
		#pad with whitespace for formatting issues in .pdb file
		#	val = join([repeat(" ",length(key)-length(val)),val]);
		#if length(key) < length(val)
		#elseif length(key) < length(val)
		#	key = join([repeat(" ",length(val)-length(key)),key]);
		#end
			outline = replace(line,key,val);
		#else
		#	outline = line;
		end
		#outline = split(outline);
		#@printf(o,"%4s%7s %5s %7s  " 
		#println(outline);
		write(o,outline);
	else
		write(o,line);
	end
end
println("Go edit formatting on testDFAImod.pdb.  Still working on this stupid bug.");
#lines = readlines(f);
-#
#run martinize to get the SC correct
run(`python martinize_hacked.py -f testDFAImod2.pdb -x CG_DFAI.pdb -n CG.ndx -nmap mapping.ndx -nt -ss ss.dat -o CG_DFAI.top -ff martini22`);
#click together central aromatics with the SCs
