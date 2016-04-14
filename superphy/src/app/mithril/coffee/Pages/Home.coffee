# coffeelint: disable=max_line_length

class Home extends Page
    Routes.add('/', @)
    Routes.add('/home', @)
    @controller: (args) ->
        args = args || {}
        @debug = args.debugging || true
        return @
    @view: (ctrl) ->
        super(
            m(".container", {id:'Home.get()-beta'},
                m("input[type=button]"
                    onclick: ->
                        localStorage.clear()
                        console.log("Cache Cleared!")
                    value: "Clear Cache"
                ) if ctrl.debug is true
                m(".row", m(".well center-block", m("p.text-center"
                    m("span.text-info beta-release", "Beta Release")
                    " Some features are still under development and may not be fully functional."
                )))
                m(".row superphy-image"
                    m("img", {src:'images/superphy_logo_with_title.png'}),
                    m("p.superphy-image"
                        'NEXT-LEVEL PHYLOGENETIC AND
                        EPIDEMIOLOGICAL ANALYSIS OF PATHOGENS'
                    )
                )
                m(".row well"
                    m("p",'A user-friendly, integrated platform for the predictive genomic analyses of '
                        m("em", 'Escherichia coli')
                    )
                    m("p", 'Provides near real-time analyses of thousands of genomes using novel computational approaches.')
                    m("p", 'Generates results that are understandable and useful to a wide community of users.')
                )
                m(".row"
                    m("button.btn btn-danger btn-lg center-block", "INTRODUCTION")
                )
                m(".row text-center"
                    m("p", 'For more information or to discuss possible collaboration, please contact:'
                        m("p", m("ul", m("li", 'Dr. Vic Gannon: vic.gannon@phac-aspc.gc.ca')))
                    )
                )
            )
        )