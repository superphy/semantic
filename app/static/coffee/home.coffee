class Home
    controller: ->

    #By default Coffeescript only returns the last item.
    #Wrap as an array to return all elements
    view: ->
        [
            header.view()
            m("div", {class:'container', id:'home-beta'}, [
                m("div", {class:'row'}, [
                    m("div", {class:'well center-block'}, [

                        m("p", {class:'text-center'}, [
                            m("span", {class:'text-info beta-release'}, "Beta Release")
                            ," Some features are still under development and may not be fully functional."
                        ])
                     ])
                ])
                m("div", {class:'row superphy-image'}, [
                    m("img", {src:'images/superphy_logo_with_title.png'}),
                    m("p", {class:'superphy-image'}, 'NEXT-LEVEL PHYLOGENETIC AND EPIDEMIOLOGICAL ANALYSIS OF PATHOGENS')
                ])
                m("div", {class:'row well'}, [
                    m("p", 'A user-friendly, integrated platform for the predictive genomic analyses of ',[
                        m("em", 'Escherichia coli')
                    ]),
                    m("p", 'Provides near real-time analyses of thousands of genomes using novel computational approaches.'),
                    m("p", 'Generates results that are understandable and useful to a wide community of users.')
                ])
                m("div", {class:'row'}, [
                    m("button", {class:'btn btn-danger btn-lg center-block'}, "INTRODUCTION")
                ])
                m("div", {class:'row text-center'}, [
                    m("p", 'For more information or to discuss possible collaboration, please contact:', [
                        m("p", [
                            m("ul", [
                                m("li", 'Dr. Vic Gannon: vic.gannon@phac-aspc.gc.ca')
                            ])
                        ])
                    ])
                ])
            ])
        ]

home = new Home()
m.mount(document["body"], home);