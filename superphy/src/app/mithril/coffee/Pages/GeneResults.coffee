# coffeelint: disable=max_line_length

class GeneResults extends Page

    Routes.add('/results', @)

    @controller: (args) ->
        @results = GeneSearchModel
        return @

    @view: (ctrl, args) ->
        ## Temp data
        testdata = {
            genomes: ["JHNV00000000", "ANVW00000000"]
            vfs: ["saa", "papC", "ompA", "hylA", "QnrS7"]
            amrs: ["QnrS7"]
        }
        data = {
            vfresults: GeneSearchModel.vfresults
            amrresults: GeneSearchModel.amrresults
        }
        return super(
            m.component(GeneResultsPanel, {
                header: m.component(ContentHeader, {
                    title: "Virulence Factor Results"
                })
                matrix: m.component(Matrix, {
                    matrixview: new MatrixView()
                    genomes: testdata.genomes
                    genes: testdata.vfs
                    elID: "genome_matrix1"
                })
                histogram: m.component(Histogram, {
                    title: "Hello World"
                })
            })
            m.component(GeneResultsPanel, {
                header: m.component(ContentHeader, {
                    title: "AMR Results"
                })
                matrix: m.component(Matrix, {
                    matrixview: new MatrixView()
                    results: data.amrresults
                    elID: "genome_matrix2"
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
        @getResults= () ->
            response = m.prop(null)
            if response() is null
                response =
                    m.request(
                        method: "POST",
                        url: "http://#{location.hostname}:5000/data/genesearchresults",
                        data: {
                            genome: args.genomes #["JHNV00000000", "ANVW00000000"]
                            genes: args.genes ## temp for testing
                        }
                        datatype: "json",
                        type: (response) ->
                            console.log("Responsee:", response)
                            return response
                    )

            return response

        @results = @getResults()
        return @

    view: (ctrl, args) ->
        m 'div', {id: 'vf_result_matrix'},
            m '.superphy-matrix', {id: args.elID, config: args.matrixview.init(ctrl.results())}


Histogram =
    view: (ctrl, args) ->
        m("h1", "Hello World")
