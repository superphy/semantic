# coffeelint: disable=max_line_length

class GeneResults extends Page

    Routes.add('/results', @)

    @controller: (args) ->
        @results = GeneSearchModel
        return @

    @view: (ctrl, args) ->
        ## Temp data for testing
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
            m '.', {id: 'wrapper'},
                m.component(Sidebar)
                m 'div', {id: "page-content-wrapper"},
                    m '.page-content inset',
                        m '.container-fluid',
                            m.component(ContentHeader, {
                                title: "Virulence Factor and AMR Results"
                            }),
                            m 'div', {id: 'results'},
                                [m.component(GeneResultsPanel, {
                                    id: "vf_results"
                                    type: "Virulence Factor"
                                    numSelected: GeneSearchModel.selectedVF.length
                                    matrix: m.component(Matrix, {
                                        matrixview: new MatrixView()
                                        results: data.vfresults
                                        parentEl: "vf_result_matrix"
                                        elID: "genome_matrix1"
                                    })
                                    histogram: m.component(Histogram, {
                                        histview: new HistogramView()
                                        results: data.vfresults
                                        parentEl: "vf_result_histogram"
                                        elID: "matrix_ticker1"
                                    })
                                }),
                                m.component(GeneResultsPanel, {
                                    id: "amr_results"
                                    type: "Antimicrobrial Resistance"
                                    numSelected: GeneSearchModel.selectedAMR.length
                                    matrix: m.component(Matrix, {
                                        matrixview: new MatrixView()
                                        results: data.amrresults
                                        parentEl: "amr_result_matrix"
                                        elID: "genome_matrix2"
                                    })
                                    histogram: m.component(Histogram, {
                                        histview: new HistogramView()
                                        results: data.amrresults
                                        parentEl: "amr_result_histogram"
                                        elID: "matrix_ticker2"
                                    })
                                })]

        )

GeneResultsPanel =
    view: (ctrl, args) ->
        if args.numSelected > 0 \
            then m('div', {id: args.id}, [
                    m('hr'),
                    m('h4', "Detected " + args.type + " Alleles")
                    args.matrix
                    args.histogram
                 ])
        else
            m('div', '')

ContentHeader =
    view: (ctrl, args) ->
        m '.row',
            m '.col-xs-8',
                m '.content-header',
                    m 'h1', args.title
            m '.col-xs-4',
                m 'button', {id: "intro-button", \
                             class: "btn btn-danger-lg", \
                             type: "button"},
                            "INTRODUCTION"

Matrix =
    controller: (args) ->
        @results = args.results()
        return @

    view: (ctrl, args) ->
        m 'div', {id: args.parentEl},
            m '.superphy-matrix', {id: args.elID, config: args.matrixview.init(ctrl.results, args.parentEl, args.elID)}


Histogram =
    controller: (args) ->
        @results = args.results()
        return @
    view: (ctrl, args) ->
        m 'div', {class: "row histogram-row"},
            m '.col-md-4 histogram-description'
                m 'span', "Blurb goes here"
            m 'div', {class: "col-md-8 histogram-container", id: args.parentEl},
                m '.superphy-histogram', {id: args.elID, config: args.histview.init(ctrl.results, args.parentEl, args.elID)}
