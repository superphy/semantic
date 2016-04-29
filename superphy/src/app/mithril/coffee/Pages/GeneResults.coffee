# coffeelint: disable=max_line_length

###
CLASS GeneResults

Page Component for the gene results page.
###
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
                m '.', {id: "page-content-wrapper"},
                    m '.page-content inset',
                        m '.container-fluid',
                            m.component(ContentHeader, {
                                title: "Virulence Factor and AMR Results"
                            }),
                            m '.', {id: 'results'},
                                ## Virulence factors
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
                                ## AMR genes
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


###
COMPONENT GeneResultsPanel

Wrapper component for results of a set of genes.

Args:
    id: ID name for div element
    type: type of gene
    numSelected: number of selected genes
    matrix: MatrixView component
    histogram: Histogram component
###
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


###
COMPONENT ContentHeader

Intro information for gene results

Args:
    title: Title of the page
###
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


###
COMPONENT Matrix

Initializer component for the Matrix

Args:
    matrixview: MatrixView component
    results: results object
    parentEl: parent element's ID
    elID: matrix's element ID
###
Matrix =
    controller: (args) ->
        @results = args.results()
        return @

    view: (ctrl, args) ->
        m 'div', {id: args.parentEl},
            m '.superphy-matrix', {id: args.elID, config: args.matrixview.init(ctrl.results, args.parentEl, args.elID)}


###
COMPONENT Histogram

Initializer component for the Histogram

Args:
    histview: Histrogram component
    results: results object
    parentEl: parent element's ID
    elID: histogram's element ID
###
Histogram =
    controller: (args) ->
        @results = args.results()
        return @
    view: (ctrl, args) ->
        m '.row histogram-row',
            m '.col-md-4 histogram-description'
                m 'span', "Blurb goes here"
            m '.col-md-8 histogram-container', {id: args.parentEl},
                m '.superphy-histogram', {id: args.elID, config: args.histview.init(ctrl.results, args.parentEl, args.elID)}
