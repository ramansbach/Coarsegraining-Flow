#A julia file to organize coarsegraining of the DFAX series
#rename all residues of .pdb file
using MATLAB
module CGWorkflow
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

function renameRes(inname,outname,resdictname,cadictname)
#renames the residues of the atomistic pdb file to match the names in the martinize script
resDict = parseDict(resdictname);
cadict = parseDict(cadictname);
println(resDict);
f = open(inname);
o = open(outname,"w");
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

function rerenameRes(CGtoAAdict,fin,fout,offset,wresname="W",v=false)
#function to take a CG gro file and rename it with atomistic residue names and original Martini beads in preparation for fine-graining
inf = open(fin);
outf = open(fout,"w");

#first of all, need to rename residues, which can be done with a dictionary that includes the number in front of the residue name
resDict = parseDict(CGtoAAdict);
rescount = 1;
for line in eachline(inf)
	spline = split(line);
	if haskey(resDict,spline[1])
		#spline[1] = resDict[spline[1]];
		if contains(spline[1],"COP") || contains(spline[1],"BEN")
			if rescount == 3
				spline[2] = "SC2";
				rescount = 1;
			else
				rescount+=1;
			end
		end
		spline[1] = resDict[spline[1]];
		num = match(r"\d+",spline[1]);
		println(spline);
		res = match(r"[a-zA-Z]+.",spline[1]);
		if v
		@printf(outf,"%5d%-4s%6s%5d % 7.3f % 7.3f % 7.3f % 7.4f % 7.4f % 7.4f\n",int(num.match),res.match,spline[2],int(spline[3]),float(spline[4]),float(spline[5]),float(spline[6]),float(spline[7]),float(spline[8]),float(spline[9])); 	
		else
		@printf(outf,"%5d%-4s%6s%5d % 7.3f % 7.3f % 7.3f\n",int(num.match),res.match,spline[2],int(spline[3]),float(spline[4]),float(spline[5]),float(spline[6])); 	
		end
	elseif contains(line,wresname)
		num = match(r"\d+",spline[1]);
		res = match(r"[a-zA-Z]+",spline[1]);
		if (length(spline) < 9 && v) || (length(spline) < 6)
			l1 = int(num.match) - offset;
			l2 = res.match;
			l3 = wresname;	
			wmatch = match(r"\d+",spline[2]);
			l4 = int(wmatch.match);
			l5 = float(spline[3]);
			l6 = float(spline[4]);
			l7 = float(spline[5]);
			if v
			l8 = float(spline[6]);
			l9 = float(spline[7]);
			l10 = float(spline[8]);
			end
		else
			l1 = int(num.match) - offset;
			l2 = res.match;
			l3 = spline[2];
			l4 = int(spline[3]);
			l5 = float(spline[4]);
			l6 = float(spline[5]);
			l7 = float(spline[6]);
			if v
			l8 = float(spline[7]);
			l9 = float(spline[8]);
			l10 = float(spline[9]);
			end
		end
		if v
		@printf(outf,"%5d%-4s%6s%5d % 7.3f % 7.3f % 7.3f % 7.4f % 7.4f % 7.4f\n",l1,l2,l3,l4,l5,l6,l7,l8,l9,l10); 	
		else
		@printf(outf,"%5d%-4s%6s%5d % 7.3f % 7.3f % 7.3f\n",l1,l2,l3,l4,l5,l6,l7); 	
	
		end

	else
		write(outf,line);
	end
end
#then rename some of the BB beads to SC ones
#then renumber the Ws (might be labeled W or PW or something) by offset 
#probably need to @printf out to get the formatting right
close(inf);
close(outf);
end
#renameRes();
#lines = readlines(f);


#run martinize to get the SC correct
#run(`bash testmart.sh`);
#run(`python martinize_hacked.py -f testDFAImod.pdb -x CG_DFAI.pdb -n CG.ndx -nmap mapping.ndx -nt -ss ss.dat -o CG_DFAI.top -ff martini22`);

#println("Currently, need to run martinize_hacked.py by hand to get this thing to work.");


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
	return bondlist;
end

function readallX(name,regex,angs)
#get all angles or dihedrals in list form from an .itp file
	f = open(name);
	bflag = false;
	bondlist = Tuple[];
	for line in eachline(f)
		 if ismatch(regex,line)
			spline = split(line);
			if angs
				bondlist = cat(1,bondlist,(int(spline[1]),int(spline[2]),int(spline[3])));
			else
				bondlist = cat(1,bondlist,(int(spline[1]),int(spline[2]),int(spline[3]),int(spline[4])));
			end
		 end
	end
	return bondlist;
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
#function th	at takes a dictionary of backbone bonds, and a dictionary of constraints, and returns an ordered set of backbone angles
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
			return reshape(angs,3,int(length(angs)/3))';
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
						nextbond = [nextind,a];
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

function getBBDihs(bondlist,clist,ind1)
#function that acts like getBBAngs, except it works on dihedrals instead
	currind = ind1;
	dihs = Int[];
	while true
		if haskey(bondlist,currind)
			bond = [currind,bondlist[currind]];
		elseif haskey(clist,currind)
			if length(clist[currind]) == 1
				bond = [currind,clist[currind]];
			else
				for d in clist[currind]
					if haskey(bondlist,d)
						bond = [currind,d];
					end
				end
			end
		else
			return reshape(dihs,4,int(length(dihs)/4))';
		end
		nextind = bond[2];
		if haskey(bondlist,nextind)
			nextbond = [nextind,bondlist[nextind]];
		elseif haskey(clist,nextind)
			if length(clist[nextind]) == 1
				nextbond = [nextind,clist[nextind]];
			else
				for d in clist[nextind]
					if haskey(bondlist,d)
						nextbond = [nextind,d];
					end
				end
			end
		else
			return reshape(dihs,4,int(length(dihs)/4))';
		end

		thdind = nextbond[2];
		if haskey(bondlist,thdind)
			thdbond = [thdind,bondlist[thdind]];
			println(bond);
			println(nextbond);
			println(thdbond);
			println("\n");
		elseif haskey(clist,thdind)
			if length(clist[thdind]) == 1
				thdind = [thdind,clist[thdind]];
			else
				for d in clist[thdind]
					if haskey(bondlist,d)
						thdbond = [thdind,d];
					end
				end
			end
		else
			return reshape(dihs,4,int(length(dihs)/4))';
		end
		dih = append!(bond,[nextbond[2]]);
		dih = append!(dih,[thdbond[2]]);
		dihs = append!(dihs,dih);
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
		remAngs[i] = (remAngs[i][1]+find,remAngs[i][2]+find,remAngs[i][3]+find);i
	end
	println(remAngs);
	bbangles = getBBAngs(bondlist,clist,ind1);
	bbdihs = getBBDihs(bondlist,clist,ind1);
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
				for i=1:size(bbangles)[1]
				
					@printf(outf,"%5d %5d %5d %6d %6d %4d ; \n",bbangles[i,1],bbangles[i,2],bbangles[i,3],2,127,20);
				end
				write(outf,line);
			elseif line=="; Backbone dihedrals\n"
				#write out the computed backbone dihedral
				write(outf,line);
				for i = 1:size(bbdihs)[1]
					@printf(outf,"%5d %5d %5d %5d %6s %6s %4s ; dummy vals\n",bbdihs[i,1],bbdihs[i,2],bbdihs[i,3],bbdihs[i,4],"X","X","X"); 
				end
				#write(outf,line);
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

function copBB(inname,outname)
#function to replace the SC2 designation in the COP residues with a BB designation because it's actually part of the backbone
	inf = open(inname);
	outf = open(outname,"w");
	for line in eachline(inf)
		if contains(line,"COP") && contains(line,"SC2")
			line = replace(line,"SC2"," BB");
		end
		write(outf,line);
	end
	close(inf);
	close(outf);
end


function input(prompt::String="")
	print(prompt)
	text = chomp(readline())
	return text
end

function getBondNumbers(filename)
#looks through an xml file to see how many bonds, angles, and dihedrals we have
	f = open(filename);
	b = 0;
	a = 0;
	d = 0;
	for line in eachline(f)
		if ismatch(r"<name>b\d+</name>",line)
			b+=1;
		elseif ismatch(r"<name>a\d+</name>",line)
			a+=1;
		elseif ismatch(r"<name>d\d+</name>",line)
			d+=1;
		end	
	end
	close(f);
	return b,a,d;
end

function inListRev(item,list)
#checks for things that are symmetric in a list--so returns the index if either the forward or backward is found
ritem = reverse(item);
if in(item,list)
	return find(list->list==item,list)[1];
elseif in(ritem,list)
	return find(list->list==ritem,list)[1];
else
	return -1;
end

end

function symmetryList(symDict)
#goes through a list of bonds, angles and dihedrals in the form (1,2) (3,4) per line (may be more than two) and returns three 
#arrays of arrays of arrays-> eg bondarray:  [ [ [1,2],[3,4] ] ] means there exists one set of symmetries including (1,2) and (3,4)
	f = open(symDict);
	bondsyms = Array[];
	angsyms = Array[];
	dihsyms = Array[];
	for line in eachline(f)
		spline = split(line);
		#println(spline);
		sym = Array[];
		for item in spline
			bond = Int[];
			spitem = split(item,['(',',',')'])[2:end-1];
			for spit in spitem
				bond = push!(bond,int(spit));
			end
			sym = push!(sym,bond);	
		end
		
		if length(sym[1]) == 2
			bondsyms = push!(bondsyms,sym);
			#println(bondsyms);
		elseif length(sym[1]) == 3
			angsyms = push!(angsyms,sym);
		elseif length(sym[1]) == 4
			dihsyms = push!(dihsyms,sym);
		end
	end
	close(f);
	return bondsyms,angsyms,dihsyms;
end


function allSymList(symDict,bbBeadOrder)
	#returns lists of all symmetries, including both backbone and non-backbone bonds, angles, dihedrals
	(bondsyms,angsyms,dihsyms) = symmetryList(symDict);
	for i = 1:int(length(bbBeadOrder)/2)-1
		bond = Array[];
		bond = push!(bond,[bbBeadOrder[i],bbBeadOrder[i+1]]);
		bond = push!(bond,[bbBeadOrder[end-i],bbBeadOrder[end-i+1]]);
		bondsyms = push!(bondsyms,bond);
		ang = Array[];
		ang = push!(ang,[bbBeadOrder[i],bbBeadOrder[i+1],bbBeadOrder[i+2]])
		ang = push!(ang,[bbBeadOrder[end-i-1],bbBeadOrder[end-i],bbBeadOrder[end-i+1]]);
		angsyms = push!(angsyms,ang);
		if i < int(length(bbBeadOrder)/2)-1
		#don't need the last one
		dih = Array[];
		dih = push!(dih,[bbBeadOrder[i],bbBeadOrder[i+1],bbBeadOrder[i+2],bbBeadOrder[i+3]]);
		dih = push!(dih,[bbBeadOrder[end-i-2],bbBeadOrder[end-i-1],bbBeadOrder[end-i],bbBeadOrder[end-i+1]]);
		dihsyms = push!(dihsyms,dih);
		end
	end
	if (length(bbBeadOrder) % 2) == 0
		#even number of backbone beads, add in last remaining bond and dihedral for the list
		midi = int(length(bbBeadOrder)/2);
		midbond = Array[];
		midbond = push!(midbond,[bbBeadOrder[midi],bbBeadOrder[midi+1]]);
		bondsyms = push!(bondsyms,midbond);
		middih = Array[];
		middih = push!(middih,[bbBeadOrder[midi-1],bbBeadOrder[midi],bbBeadOrder[midi+1],bbBeadOrder[midi+2]]);
		dihsyms = push!(dihsyms,middih);
	end
	return bondsyms,angsyms,dihsyms;
end

function writeVotcaInput(symDict,bbBeadOrder,T,outfilename,topname)
	#function that computes the list of bonds to write out to a VOTCA file, matches those to their corresponding index, and writes out a VOTCA input file to compute all (symmetric) bond data
        angreg = r"\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+\d+\s+;";
	dihreg = r"\s+\d+\s+\d+\s+\d+\s+\d+\s+[0-9a-zA-Z]+\s+[0-9a-zA-Z]+\s+[0-9a-zA-Z]+\s+;";	
	bondsyms,angsyms,dihsyms = allSymList(symDict,bbBeadOrder);
 	allbonds = readallbonds(topname);	#these bonds come in as a list of tuples, need to be careful when checking
	allangs = readallX(topname,angreg,true);     #still a list of tuples
	#println(allangs);
	alldihs = readallX(topname,dihreg,false);     #continuing to be a list of tuples
	println(alldihs);
	outf = open(outfilename,"w");
	@printf(outf,"tab set T %3d\n",T);
	@printf(outf,"tab set scale bond\n");
	sbind = 1;
	tabstrb = "tab bsX.pot";
	tabstra = "tab asX.pot";
	tabstrd = "tab dsX.pot";
	hstrb = "hist hbsX.dat";
	hstra = "hist hasX.dat";
	hstrd = "hist hdsX.dat";
	vstrb = "vals vbsX.dat";
	vstra = "vals vasX.dat";
	vstrd = "vals vdsX.dat";
	for symbond in bondsyms
		tabout = replace(tabstrb,"X",string(sbind));
		hout = replace(hstrb,"X",string(sbind));
		vout = replace(vstrb,"X",string(sbind));
		for bond in symbond
			bind = inListRev((bond[1],bond[2]),allbonds);
			if bind == -1
				println("Error.  Bond ",bond," in symmetric bond list not found in full bond list.  Check bonds.");
				return;
			else
				tabout = tabout*" *:bX:*";
				tabout = replace(tabout,"X",string(bind));
				hout = hout*" *:bX:*";
				hout = replace(hout,"X",string(bind));
				vout = vout*" *:bX:*";
				vout = replace(vout,"X",string(bind));
			end
			
		end
		write(outf,tabout*"\n");
		write(outf,hout*"\n");
		write(outf,vout*"\n");
		sbind+=1;
	end
	saind = 1;
	@printf(outf,"tab set scale angle\n");
	for symang in angsyms
		tabout = replace(tabstra,"X",string(saind));
		hout = replace(hstra,"X",string(saind));
		vout = replace(vstra,"X",string(saind));
		for ang in symang
			aind = inListRev((ang[1],ang[2],ang[3]),allangs);
			if aind == -1
				println("Error.  Angle ",ang," in symmetric angle list not found in full angle list.  Check angles.");
				return;
			else
				tabout = tabout*" *:aX:*";
				tabout = replace(tabout,"X",string(aind));
				hout = hout*" *:aX:*";
				hout = replace(hout,"X",string(aind));
				vout = vout*" *:aX:*";
				vout = replace(vout,"X",string(aind));
			end
			
		end
		write(outf,tabout*"\n");
		write(outf,hout*"\n");
		write(outf,vout*"\n");
		saind+=1;
	end
	sdind = 1;
	@printf(outf,"tab set scale no\n");
	for symdih in dihsyms
		tabout = replace(tabstrd,"X",string(sdind));
		hout = replace(hstrd,"X",string(sdind));
		vout = replace(vstrd,"X",string(sdind));
		for dih in symdih
			dind = inListRev((dih[1],dih[2],dih[3],dih[4]),alldihs);
			if dind == -1
				println("Error.  Dihedral ",dih," in symmetric dihedral list not found in full dihedral list.  Check dihedrals.");
				return;
			else
				tabout = tabout*" *:dX:*";
				tabout = replace(tabout,"X",string(dind));
				hout = hout*" *:dX:*";
				hout = replace(hout,"X",string(dind));
				vout = vout*" *:dX:*";
				vout = replace(vout,"X",string(dind));
			end
			
		end
		write(outf,tabout*"\n");
		write(outf,hout*"\n");
		write(outf,vout*"\n");
		sdind+=1;
	end
	close(outf);	
end

#need to go to votca and be able to call  (cat inputfile.txt) | csg_boltzmann --top md.tpr --trj md.xtc --cg "CGDFAG.xml;water.xml"
#that means the CG.xml file needs to be correct.  make_xml.py is supposed to do this, but I don't think it's set up to fully 
#do the dihedral angles as well.  So the first step is to get a functional xml file.  Then we can run csg_boltzmann
#We'll also need to create an inputfile (it should be very similar to the one for DFAG)
#okay so make_xml still needs to be set up to do backbone dihedrals, and we need to pay attention to the bond order
#because it's different from what I expected and that may throw things off if I don't track it properly--oh it's because we changed
#the bond order when we were rebonding? Ah, yes, that's right, so that should be fine, because the .xml file will match the CG file
#Just bear it in mind.  But yeah if we can just add BB dihedrals, we'll have a solid .xml file
#We should just be able to write the dihedrals into the CG .top file and that'll do well enough?
function main(prompts=true)
	#set default values
	prepform = "y";
	pdbname = "testDFAI";
	rdictname = "resDict.txt";
	cdictname = "CAdict.txt";
	itpcorr = "y";
	ind1 = 10;
	outxml = "CG_DFAI.xml";
	AAtop = "DFAI.top";
	CGtop = "Protein_rebond_reangle_reBB.itp";
	map = "mapping.ndx";
	bbs = [1,3,7,8,10,11,13,14,15,17,18,22,20,19,23,25,26,30]; #list of back-bone beads in order, currently hardcoded for DFAI
	symdictname = "DFAIsymdict.txt";
	votcainputname = "DFAIvotcasyminput.txt";
	mdtop = "/home/rachael/cluster/scratch/ramansbach/coarsegraining/AA/DFAI/4_md/md.tpr";
	mdtrj = "/home/rachael/cluster/scratch/ramansbach/coarsegraining/AA/DFAI/4_md/md.xtc";
	T=298; #temperature to run VOTCA at
	#function to step through the process
	#currently prompts the user, need to add a version that runs from an input file
	if prompts
	println("Welcome to the DFAX parser.  This should step you through the process of
	coarsegraining a molecule of the DFAX series, as well as testing the 
	model you get, and setting up to do long MD runs.  First, we need to 
	prepare an atomistic .pdb file for running through martinize.py.  \n")
	end
	
	if prompts prepform = input("Would you like to do this now? y/n "); end
	if prepform == "y"
		pdbname = input("What is the name of the pdb file? (without extension) ");
		inname = pdbname*".pdb";
		outname = pdbname*"mod.pdb";
		renameRes(inname,outname,rdictname,cdictname);		
	end
	if prompts 
		println("Now that the pdb file is properly prepared, we can run martinize.py
		For the time being, this has to be done by hand.  Please run:\n
		python martinize_hacked.py -f testDFAImod.pdb -x CG_DFAI.pdb -n CG.ndx -nmap mapping.ndx -nt -ss ss.dat -o CG_DFAI.top -ff martini22\n
		where the inputs are testDFAImod.pdb (the input name of the pdb file + mod), and ss.dat
		which must be in the folder.  For polarizable water, specify martini22p instead.\n");
	end
	if prompts 
	println("The .pdb file is correct, but the bonds and angles are not.  They have to be corrected.");
	itpcorr = input("Would you like to do this now? y/n ");
	end
	if itpcorr == "y"
		ind1 = int(input("What number residue marks the last residue before the aromatics? "));
		reBond("Protein.itp","Protein_rebond.itp",ind1);
		r1 = int(input("What is the index of the first backbone residue? "));
		reAngle("Protein_rebond.itp","Protein_rebond_reangle.itp",r1,ind1);
		copBB("Protein_rebond_reangle.itp","Protein_rebond_reangle_reBB.itp");				
	end
	if prompts
		println("The next step is to go to VOTCA and create potentials and histograms of
			all bonds of interest.");
		xmlrun = input("First we must prepare an .xml topology file for use with VOTCA.  Would you like to do this now? y/n");
		
	end
	if xmlrun == "y"
		fnames = input("Would you like to input non-default file names? y/n");
		if fnames == "y"
			outxml = input("What is the name of the output xml file? ");
			AAtop = input("What is the name of the atomistic topology file? ");
			CGtop = input("What is the name of the coarsegrained itp file? ");
			map = input("What is the name of the mapping index file? ");
		end
		run(`python ../../make_xml.py -o $outxml -a $AAtop -c $CGtop -m $map `);
	end
	if prompts
		println("Now we must prepare the input file that tells VOTCA which potentials to output.");
		prepin = input("Would you like to do this now? y/n ");
	end
	if prepin == "y"
		#bn,an,dn = getBondNumbers(outxml);
		writeVotcaInput(symdictname,bbs,T,votcainputname,"Protein_rebond_reangle_reBB.itp");
	end
	if prompts
		println("Now that the input file is prepared, we can run VOTCA using the produced input and xml files.  This will produce a lot of data files and require knowing the locations of the trajectory and tpr files we want to analyze.");
		runvotca = input("Would you like to do this now? y/n ");
	end
	if runvotca == "y"
		giveinputfilenames = input("Would you like to give the locations and names of the input files? y/n ");
		if giveinputfilenames=="y"
			println("This functionality is not yet set-up.");
		else
		#need to go to votca and be able to call  (cat inputfile.txt) | csg_boltzmann --top md.tpr --trj md.xtc --cg "CGDFAG.xml;water.xml"
			#println("Currently this functionality cannot be run from within the program.  Please run\n
			#cat <votcainputname>) | csg_boltzmann --top <mdtop> --trj <mdtrj> --cg \"<outxml>;water.xml\"");
			#run(`(cat $votcainputname) | csg_boltzmann --top $mdtop --trj $mdtrj --cg '$outxml;water.xml'`);
			mat"runvotca(votcainputname,mdtop,mdtrj,outxml)";
		end
	end

end
end
#main()
