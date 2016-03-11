class Factors extends Page_Template
    #a list to be passed to the link generator function
    model: () => 
        @value = m.prop('')
        
        @table ?= new Table()
        @search ?= new GeneSearch()
        @vfform ?= new GeneForm('Virulence Factor')
        @amrform ?= new GeneForm('Antimicrobial Resistance')
        @sidebar ?= new Sidebar()
        return

    controller: (options) =>
        @model()
        @data ?= new GeneData()
        return
        # data = @data
        # return {
        #     pages: pages,
        #     rotate: ->
        #         pages().push(pages().shift)
        # }
    view: (ctrl) =>
        return [
            super()
            m('div', {id: 'wrapper'}, [
                @sidebar.view()
                m('div', {id: 'page-content-wrapper'}, [
                    m('div', {id: 'page-content -inset'}, [
                        m('div', {class: 'container-fluid'}, [
                            m('div', {class: 'row'}, [
                                m('div', {class: 'col-xs-8'}, [
                                    m('div', {class: 'intro'}, [
                                        m('p', 'BSearch for the presence or absence of virulence factor genes or antimicrobial resistance 
                                                genes in genomes of interest. Detailed information on individual virulence factor or 
                                                antimicrobial resistance genes can be retrieved by clicking on the individual genes.')
                                    ])

                                    m('div', {class: 'search'}, [
                                        m('ul', {id: 'gene-search-tabs', class: 'nav nav-tabs nav-justified'}, [
                                            m('li', {class: 'active'}, m("a", {href:"#gene-search-querygenes", id: "select_genes", config:m.route}, "Select Genes")),
                                            m('li', m("a", {href:"#gene-search-genomes", id: "select_genomes", config:m.route}, "Select Genomes")),
                                            m('li', m("a", {href:"#gene-search-submit", id: "submit_query", config:m.route}, "Submit Query"))
                                        ])
                                    ])

                                    m('div', {class: 'tab-content'}, [
                                        m('div', {class: 'tab-pane active', id: 'gene-search-querygenes'}, [
                                            m('div', {class: 'panel-group genes-search', id: 'accordian'}, [
                                                @vfform.view()
                                                @amrform.view()
                                            ])
                                        ])
                                    ])
                                    @search.view(@data)
                                    # @table.view(@data)
                                ])
                            ])
                        ])
                    ])
                ])
            ])
        ]