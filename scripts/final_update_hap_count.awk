# using the unique haplotypes to get the modified haplotype counts for each population (thus accounting missing haplotypes)

NR==FNR{
    uniq[$1]++;
    next
}

{
    if(uniq[$1]==1)
    {
	A[$1]=$2
	uniq[$1]=0
    }
}

END{
    for(item in uniq)
    {
	if(uniq[item]==1)
	{
	    A[item]=0
	}
    }
    for(hap in A)
    {
	printf ("%s:%d\t",hap, A[hap])
    }
    printf("\n")
}
