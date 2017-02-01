from rdflib import BNode, Literal, Graph

# working with Serotype, Antimicrobial Resistance, & Virulence Factor data structures
def parse_serotype(graph, serotyper_dict, uriIsolate):
    if 'O type' in serotyper_dict:
        graph.add((uriIsolate, gu('ge:0001076'),
                   Literal(serotyper_dict['O type'])))
    if 'H type' in serotyper_dict:
        graph.add((uriIsolate, gu('ge:0001077'),
                   Literal(serotyper_dict['H type'])))
    if 'K type' in serotyper_dict:
        graph.add((uriIsolate, gu('ge:0001684'),
                   Literal(serotyper_dict['K type'])))

    return graph


def parse_gene_dict(graph, gene_dict, uriGenome):
    '''
    My intention is to eventually use ECTyper for all of the calls it was meant for.
    Just need to update ECTyper dict format to ref. AMR / VF by contig. as opposed to genome directly.

    These are the common gene related triples to both AMR / VF.
    Note: we are working from uriGenome and assume that the calling functions (
    generate_amr() and generate_vf() are doing the transformations to the
    gene_dict.keys so that they are contig ids (as they differ in return value
    between VF & AMR from ECTyper)
    )

    TODO: offshore rgi calls to ectyper and make it return a dict in the format we need
    -currently, we'll handle ORF_ID to contig id transform in generate_amr()

    Args:
    graph(rdflib.Graph): the running graph with all our triples
    gene_dict({{}}): a dictionary of genes with a assoc info
        ex. {'Some_Contig_ID':[{'START','STOP','ORIENTATION','GENE_NAME'}]}
    uriGenome(rdflib.URIRef): the base uri of the genome
        ex. :4eb02f5676bc808f86c0f014bbce15775adf06ba

    TODO: merge common components with generate_amr()
    '''

    for contig_id in gene_dict.keys():
        for gene_record in gene_dict[contig_id]:

            # recreating the contig uri
            uriContig = gu(':' +
                           contig_id)  # now at contig uri

            # after this point we switch perspective to the gene and build down to
            # relink the gene with the contig

            bnode_start = BNode()
            bnode_end = BNode()

            # some gene names, esp those which are effectively a description,
            # have spaces
            gene_name = gene_record['GENE_NAME'].replace(' ', '_')

            graph.add((gu(':' + gene_name), gu('faldo:Begin'), bnode_start))
            graph.add((gu(':' + gene_name), gu('faldo:End'), bnode_end))

            # this is a special case for amr results
            if 'CUT_OFF' in gene_dict.keys():
                graph.add((bnode_start, gu('dc:Description'),
                           Literal(amr_results['CUT_OFF'][i])))
                graph.add((bnode_end, gu('dc:Description'),
                           Literal(amr_results['CUT_OFF'][i])))

            graph.add((bnode_start, gu('rdf:type'), gu('faldo:Position')))
            graph.add((bnode_start, gu('rdf:type'), gu('faldo:ExactPosition')))
            graph.add((bnode_end, gu('rdf:type'), gu('faldo:Position')))
            graph.add((bnode_end, gu('rdf:type'), gu('faldo:ExactPosition')))

            if gene_record['ORIENTATION'] is '+':
                graph.add((bnode_start, gu('rdf:type'), gu(
                    'faldo:ForwardStrandPosition')))
                graph.add((bnode_end, gu('rdf:type'), gu(
                    'faldo:ForwardStrandPosition')))
            else:
                graph.add((bnode_start, gu('rdf:type'), gu(
                    'faldo:ReverseStrandPosition')))
                graph.add((bnode_end, gu('rdf:type'), gu(
                    'faldo:ReverseStrandPosition')))

            graph.add((bnode_start, gu('faldo:Position'),
                       Literal(gene_record['START'])))
            graph.add((bnode_start, gu('faldo:Reference'), uriContig))

            graph.add((bnode_end, gu('faldo:Position'),
                       Literal(gene_record['STOP'])))
            graph.add((bnode_end, gu('faldo:Reference'), uriContig))

            ####

    return graph
