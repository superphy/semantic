# coffeelint: disable=max_line_length

class GeneResults extends Page
    @controller: (args) ->
        return @

    @view: (ctrl, args) ->
        ## Temp data
        data = {
            genomes: ["JHNV00000000", "ANVW00000000"]
            genes: ["saa", "papC", "ompA", "hylA"]
        }
        return super(
            m.component(GeneResultsPanel, {
                header: m.component(ContentHeader, {
                    title: "Virulence & AMR Results"
                })
                matrix: m.component(Matrix, {
                    matrixview: new MatrixView()
                    genomes: data.genomes
                    genes: data.genes
                    elID: "genome_matrix1"
                })
                histogram: m.component(Histogram, {
                    title: "Hello World"
                })
            })

        )

GeneResultsPanel =
    view: (ctrl, args) ->
        m 'div', {id: "page-content-wrapper"},
            m '.page-content inset',
                m '.container-fluid',
                    args.header,
                    m 'div', {id: 'results'},
                        m 'div', {id: 'vf_results'},
                            m('hr'),
                            m('h4', "Detected Virulence Factor Alleles"),
                            args.matrix,
                            args.histogram

ContentHeader =
    view: (ctrl, args) ->
        m '.row',
            m 'col-xs-8',
                m '.content-header',
                    m 'h1', args.title
            m 'col-xs-4'
                m 'button', {id: "intro-button", \
                             class: "btn btn-danger-lg", \
                             type: "button"},
                            "INTRODUCTION"

Matrix =
    controller: (args) ->
        self = @
        @genomeDict = {}
        for genome in args.genomes
            @genomeDict[genome] = {}
            for gene in args.genes
                @genomeDict[genome][gene] = 0

        request = (selectedGenome, json={}) ->
            parseResponse=(response)->
                for binding in response.results.bindings
                    gene_name = binding["Gene_Name"]['value']
                    self.genomeDict[selectedGenome][gene_name] += 1
                console.log("Response:", response)

            #ep = getEndpoint("data/genesearchresults")

            m.request(
                method: "POST",
                url: "http://#{location.hostname}:5000/data/genesearchresults",
                data: {
                    genome: selectedGenome
                    genes: ["saa"] ## temp for testing
                }
                datatype: "json",
                type: parseResponse
            )

        for g in args.genomes
            request(g)

        return @

    view: (ctrl, args) ->
        m 'div', {id: 'vf_result_matrix'},
            m '.superphy-matrix', {id: args.elID, config: args.matrixview.init(ctrl.genomeDict)}


Histogram =
    view: (ctrl, args) ->
        m("h1", "Hello World")