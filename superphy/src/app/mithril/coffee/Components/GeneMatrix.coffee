class GeneMatrix
    constructor: (@genomes, @genes, @results) ->
        ## This might change, but for now, @genomes is a list of genome accessions,
        ## @genes is a list of gene names, and
        ## @results are the results returned from the query.

        ## Set-up for d3
        @cellWidth = 20
        @margin = {top: 150, right: 0, bottom: 0, left:250}
        @height = @genomes.length 

    computeMatrix: (genomes, genes, alleles) ->

