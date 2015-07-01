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


function reBond(inname,outname,ind1)
#pull in .itp file, remove all incorrect bonds (with indices remapped)
#also add new correct bonds (indices still remapped)
#currently only works for DFAX series--may want to expand for aromatics in general?
	remBonds = [(1,2),(5,6),(10,11),(1,4),(4,5),(5,8),(8,9),(9,10),(10,13),(2,3),(2,4),(6,7),(6,8),(11,12),(11,13)];
	addBB = reshape([3,4,4,5,7,8,8,12,10,9,9,13],2,6)';
	addC = reshape([1,2,1,3,2,3,5,7,5,6,6,7,12,10,12,11,11,10],2,9)';
	for i = 1:length(remBonds)
		remBonds[i] = (remBonds[i][1]+ind1,remBonds[i][2]+ind1);
	end
	println(remBonds);
	addBB+=ind1;
	addC+=ind1;
	inf = open(inname);
	outf = open(outname,"w");
	for line in eachline(inf)
		if line=="; Backbone bonds\n"
			write(outf,line);
			for i = 1:size(addBB)[1]
			
				@printf(outf,"%5d %5d %6d   %5.5f %5d ; \n",addBB[i,1],addBB[i,2],1,0.35,1250);
				
			end
		elseif line=="[ constraints ]\n"
			write(outf,line);
			#println(addC);
			for i = 1:size(addC)[1]
				#println(i)
				#println((addC[i][1],addC[i][2]));
				@printf(outf,"%5d %5d %6d   %5.5f ; \n",addC[i,1],addC[i,2],1,0.27);
			end
		
		else
			
			if ismatch(r"\s+\d+\s+\d+\s+\d+\s+\d+\.\d+\s+;",line) || ismatch(r"\s+\d+\s+\d+\s+\d+\s+\d+\.\d+\s+\d+\s+;",line)
				spline = split(line)
				bond = (int(spline[1]),int(spline[2]));
				#println("bond is: ");
				#println(bond);
				#println("searching in: ");
				#println(remBonds);
				#println("checking",line);
				#println("bond is: ",bond);
				#println("remBonds is:", remBonds);
				remBond = in(bond,remBonds);
				if !(remBond)
					write(outf,line);
				else
					println("not writing out",line);
				end
			else
				write(outf,line);
			end
	
		end

	end
close(inf);
close(outf);
end
#fix bonds/angles/dihedrals in the .itp file
