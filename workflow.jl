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

function readallbonds(name)
#get all bonds in list form from an .itp file
	f = open(name);
	bflag = false;
	bondlist = Tuple[];
	for line in eachline(f)
		 if ismatch(r"\s+\d+\s+\d+\s+\d+\s+\d+\.\d+\s+;",line) || ismatch(r"\s+\d+\s+\d+\s+\d+\s+\d+\.\d+\s+\d+\s+;",line)
			spline = split(line);
			bondlist = cat(1,bondlist,(int(spline[1]),int(spline[2])));
		 end
	end
	return bondlist;i
end

function readbonds(name,readstart,readend)
#get the bonds from an .itp file and return them as a dictionary
	f = open(name);
	bbflag = false;
	bondlist = Dict();
	for line in eachline(f)
		if line == readstart
			bbflag = true;
		elseif line == readend
			close(f);
			return bondlist;
		end
		if bbflag && line!= readstart && line!="\n"
			spline = split(line);
			#println(spline);
			if haskey(bondlist,int(spline[1]))
				bondlist[int(spline[1])] = append!([bondlist[int(spline[1])]],[int(spline[2])])
			else
				bondlist[int(spline[1])] = int(spline[2]);
			end	
		else
		end
	end
	close(f);
	return bondlist;
end


function readconstraints(name)
#get the constraints from an .itp file and return them as a list of tuples
	f = open(name);
	cflag = false;
	clist = Dict();
	for line in eachline(f)
		if line == "[ constraints ]\n"
			cflag = true;
		elseif line == "[ angles ]\n"
			close(f);
			return clist;
		end
		if cflag && line!="\n" && line!= "[ constraints ]\n"
			spline = split(line);
			if haskey(clist,int(spline[1]))
				clist[int(spline[1])] = append!([clist[int(spline[1])]],[int(spline[2])]);
			else
				clist[int(spline[1])]=int(spline[2]);
			end	
		else
		end
	end
	close(f);
	return clist;

end


function getBBAngs(bondlist,clist,ind1)
#function that takes a list of backbone bond tuples, and a list of constraints, and returns an ordered set of backbone angles
	currind = ind1;
	angs = Int[];
	while true
		if haskey(bondlist,currind)
			bond = [currind,bondlist[currind]];
		elseif haskey(clist,currind)
			if length(clist[currind])==1
				bond = [currind,clist[currind]];
			else
				for a in clist[currind]
					if haskey(bondlist,a)
						bond = [currind,a];
					end
				end
			end
		else
			#println(angs);
			return reshape(angs,3,int(length(angs)/3))';i
		end
		nextind = bond[2];
		if haskey(bondlist,nextind)
			nextbond = [nextind,bondlist[nextind]];
		elseif haskey(clist,nextind)
			if length(clist[nextind])==1
				nextbond = [nextind,clist[nextind]];
			else
				for a in clist[nextind]
					if haskey(bondlist,a)
						nextbond = [currind,a];
					end
				end
			end
		else
			#println(angs);
			return reshape(angs,3,int(length(angs)/3))';
		end
		ang = append!(bond,[nextbond[2]]);
		#println("ang is",ang);
		angs = append!(angs,ang);
		currind = nextind;
	end
end

function reverseDict(d)
#swaps key-value pairs to value-key pairs and returns a "reversed" dictionary
	rd = Dict();
	for key in keys(d)
		rd[d[key]] = key;
	end
	return rd;
end

function reAngle(inname,outname,ind1,find,DEBUG=false)
#redo BB angles given new BB bonds
#remove set of defined aromatic angles that are definitively wrong (hardcoded for DFAX currently)
#replace illegal bonds in angles with at least one legal bond with the appropriate legal bond
	bondlist = readbonds(inname,"; Backbone bonds\n","; Sidechain bonds\n");
	allbondlist = readallbonds(inname);
	sclist = readbonds(inname,"; Sidechain bonds\n","[ constraints ]\n");
	clist = readbonds(inname,"[ constraints ]\n","[ angles ]\n");
	remAngs = [(0,1,2),(4,5,6),(9,10,11),(1,2,3),(5,6,7),(10,11,12)];
	for i=1:length(remAngs)
		remAngs[i] = (remAngs[i][1]+find,remAngs[i][2]+find,remAngs[i][3]+find);
	end
	println(remAngs);
	bbangles = getBBAngs(bondlist,clist,ind1);
	inf = open(inname);
	outf = open(outname,"w");
        pastBBflag = false;
	
	for line in eachline(inf)
		if line=="; Backbone angles\n"
			write(outf,line);
			pastBBflag = true;
		elseif line=="; Backbone-sidechain angles\n"
			pastBBflag = false;
		end
		if !pastBBflag
			if line=="; Backbone-sidechain angles\n"
				#write out the computed backbone angles before writing out the divider
				for i=1:size(bbangles)[1]-1
				
					@printf(outf,"%4d %5d %5d %6d %6d %4d ; \n",bbangles[i,1],bbangles[i,2],bbangles[i,3],2,127,20);
				end
				write(outf,line);
			elseif ismatch(r"\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+;",line)
				#if we're past the angles we know we dont' want to write out and we're an angle or dihedral
				
				spline = split(line);
				if DEBUG println("spline is: ",spline) end
				#if it's a dihedral, it's fine, write it out
				if spline[8] == ";"
					write(outf,line);
				
				
				elseif !(in((int(spline[1]),int(spline[2]),int(spline[3])),remAngs))
					#if it's not one of the angles forbidden by the DFAX series then we can look at it more
					if DEBUG println("Not forbidden by DFAX."); end
					bond1 = (int(spline[1]),int(spline[2]));
					bond2 = (int(spline[2]),int(spline[3]));
					if (in(bond1,allbondlist) || in((bond1[2],bond1[1]),allbondlist)) && (in(bond2,allbondlist)||in((bond2[2],bond2[1]),allbondlist))
					#if both the bonds making up the angle are legal, write the angle out
						if DEBUG println("Totally legal angle"); end
						write(outf,line);
					elseif (in(bond1,allbondlist) || in((bond1[2],bond1[1]),allbondlist))
						if DEBUG println("Bond 1 legal, Bond 2 illegal."); end
					#if the first bond is legal, try to fix it
						if haskey(bondlist,bond1[2])
							if DEBUG println("Bond 1 in BB bonds, applying bond 2."); end
							#if the first bond is a backbone bond, find the matching second bond and write it out
							#Would need to add more code here if I wanted to be more careful about parameters, but for now we're not worrying about angle parameters (because they're all tabulated), so as long as the format is right, we're fine
							@printf(outf,"%4d %5d %5d %6d %6d %4d ; \n",bond1[1],bond1[2],bondlist[bond1[2]],2,127,20);
						elseif haskey(reverseDict(bondlist),bond1[2])
							if DEBUG println("Bond 1 in BB bonds, applying reversed bond 2."); end
							@printf(outf,"%4d %5d %5d %6d %6d %4d ; \n",bond1[1],bond1[2],reverseDict(bondlist)[bond1[2]],2,127,20);
						elseif haskey(sclist,bond1[2])
							if DEBUG println("Bond 1 in SC bonds, applying bond 2"); end		
							#or the first bond might be a SC bond or a constraint
							@printf(outf,"%4d %5d %5d %6d %6d %4d ; \n",bond1[1],bond1[2],sclist[bond1[2]],2,127,20);
						elseif haskey(reverseDict(sclist),bond1[2])
							if DEBUG println("Bond 1 in SC bonds, applying reversed bond 2"); end		
							#or the first bond might be a SC bond or a constraint
							@printf(outf,"%4d %5d %5d %6d %6d %4d ; \n",bond1[1],bond1[2],reverseDict(sclist)[bond1[2]],2,127,20);
						elseif haskey(clist,bond1[2])
							if DEBUG println("Bond 1 in constraints, applying bond 2"); end
							#need a try/catch here in case clist has more than one bond matching that bond
							try
								@printf(outf,"%4d %5d %5d %6d %6d %4d ; \n",bond1[1],bond1[2],clist[bond1[2]],2,127,20);
							catch
								println("Error: Possibly ambiguous angle.  Unclear secondary bond.");
							end
						elseif haskey(reverseDict(clist),bond1[2])
							if DEBUG println("Bond 1 in constraints, applying bond 2"); end
							#need a try/catch here in case clist has more than one bond matching that bond
							try
								@printf(outf,"%4d %5d %5d %6d %6d %4d ; \n",bond1[1],bond1[2],reverseDict(clist)[bond1[2]],2,127,20);
							catch
								println("Error: Possibly ambiguous angle.  Unclear secondary bond.");
							end
						else
						#if the first bond isn't anywhere we're done
						if DEBUG println("Bond 1 not found."); end
							continue;
						end
					elseif (in(bond2,allbondlist) || in((bond2[2],bond2[1]),allbondlist))
						#second bond is legal, first is illegal
						if DEBUG println("Bond 2 legal, bond 1 illegal."); end
						if haskey(bondlist,bond2[1])
							if DEBUG println("Bond 2 in BB bonds, applying reverse bond 1"); end
						#if the 2nd bond is a backbone bond, find the matching first bond and write it out
						#Would need to add more code here if I wanted to be more careful about parameters, but for now we're not worrying about angle parameters (because they're all tabulated), so as long as the format is right, we're fine
							@printf(outf,"%5d %5d %5d %6d %6d %4d ; \n",bondlist[bond2[1]],bond2[1],bond2[2],2,127,20);
						elseif haskey(sclist,bond2[1])
							if DEBUG println("Bond 2 in SC bonds,applying reverse bond 1"); end		
						#or the 2nd bond might be a SC bond or a constraint
							@printf(outf,"%5d %5d %5d %6d %6d %4d ; \n",sclist[bond2[1]],bond2[1],bond2[2],2,127,20);
						elseif haskey(clist,bond2[1])
							if DEBUG println("Bond 2 in constraints, applying reverse bond 1"); end
							#need a try/catch here in case clist has more than one bond matching that bond
							try
								@printf(outf,"%5d %5d %5d %6d %6d %4d ; \n",clist[bond2[1]],bond2[1],bond2[2],2,127,20);
							catch
								println("Error: Possibly ambiguous angle.  Unclear secondary bond.");
							end
						elseif haskey(reverseDict(bondlist),bond2[1])
							if DEBUG println("Bond 2 in BB bonds."); end
						#if the 2nd bond is a backbone bond, find the matching first bond and write it out
						#Would need to add more code here if I wanted to be more careful about parameters, but for now we're not worrying about angle parameters (because they're all tabulated), so as long as the format is right, we're fine
							@printf(outf,"%5d %5d %5d %6d %6d %4d ; \n",reverseDict(bondlist)[bond2[1]],bond2[1],bond2[2],2,127,20);
						elseif haskey(reverseDict(sclist),bond2[1])
							if DEBUG println("Bond 2 in SC bonds."); end		
						#or the 2nd bond might be a SC bond or a constraint
							@printf(outf,"%5d %5d %5d %6d %6d %4d ; \n",reverseDict(sclist)[bond2[1]],bond2[1],bond2[2],2,127,20);
						elseif haskey(reverseDict(clist),bond2[1])
							if DEBUG println("Bond 2 in constraints."); end
							#need a try/catch here in case clist has more than one bond matching that bond
							try
								@printf(outf,"%5d %5d %5d %6d %6d %4d ; \n",reverseDict(clist)[bond2[1]],bond2[1],bond2[2],2,127,20);
							catch
								println("Error: Possibly ambiguous angle.  Unclear secondary bond.");
							end
						#if the 2nd bond isn't anywhere we're done
						else
							if DEBUG println("Bond 2 not found."); end
							continue;
						end
					else
						if DEBUG println("Both bonds illegal, throw out angle."); end
						#both bonds are illegal, we may as well junk the angle
						continue;
					end
				end

				#end
			else
				println("Not a bond.");
				write(outf,line);
			end
			
		end
	end	
	close(inf);
	close(outf);
end
