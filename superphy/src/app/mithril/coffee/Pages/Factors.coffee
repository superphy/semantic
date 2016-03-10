class Factors extends Page_Template
    #a list to be passed to the link generator function
    model: () => 
        @value = m.prop('')
        @data ?= new GeneData()
        @table ?= new Table()
        @search ?= new GeneSearch()
        @vfform ?= new GeneForm('Virulence Factor')
        @amrform ?= new GeneForm('Antimicrobial Resistance')
        # return
        # @data = [
        #     {id: 1, url: "John"}
        #     {id: 2, url: "Mary"}
        #     {id: 3, url: "Jimin"}
        # ]
        return

    controller: (options) =>
        @model()
        return
        # data = @data
        # return {
        #     pages: pages,
        #     rotate: ->
        #         pages().push(pages().shift)
        # }
    view: () =>
        return [
            super()
            m('div', {class: 'intro'}, [
                m("p", 'BSearch for the presence or absence of virulence factor genes or antimicrobial resistance 
                        genes in genomes of interest. Detailed information on individual virulence factor or 
                        antimicrobial resistance genes can be retrieved by clicking on the individual genes.')])

            m('div', {class: 'search'}, [
                m('ul', {id: 'gene-search-tabs', class: 'nav nav-pills nav-justified'}, [
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
        ]
        # return[
        #     super()
        #     m("div", [
        #         # ctrl.pages().map((page) ->
        #         #     m('a', {href: page.url}, page.id)
        #         #     )
        #         m('input', {oninput: m.withAttr('value', @value), value : @value() })
        #         m('button', {onclick: => @data.search(@value())}, "Rotate links")
        #     ])
        # ]

    # # model: () =>
    # #     @value = m.prop('')
    # #     @data ?= new Data()
    # #     @table ?= new Table()
    # #     return
    # controller: () =>
    #     # for testing purposes
    #     @data = [
    #         {
    #             id: 1
    #             name: 'John'
    #         }
    #         {
    #             id: 2
    #             name: 'Mary'
    #         }
    #         {
    #             id: 3
    #             name: 'Seniqua'
    #         }
    #     ]
    #     @currentUser = m.prop(data[1])
    #     @changeUser = (id) =>
    #         console.log(id)
    #     return 
    # view: () =>
    #     return [
    #         super()
    #         m('div', {class: 'intro'}, [
    #             m("p", 'BSearch for the presence or absence of virulence factor genes or antimicrobial resistance genes in genomes of interest. Detailed information on individual virulence factor or antimicrobial resistance genes can be retrieved by clicking on the individual genes.')])
    #         # m("div", [
    #         #     m("label", "User:"),
    #         #     m.component(Select, {
    #         #         data: ctrl.data,
    #         #         value: ctrl.currentUser
    #         #         onchange: ctrl.changeUser
    #         #     })
    #         # ])
    #     ]