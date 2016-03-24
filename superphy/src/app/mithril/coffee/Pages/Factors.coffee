class Factors extends PageTemplate
    #a list to be passed to the link generator function

    model: () =>
        @value = m.prop('')
        @activeTab = m.prop('genes')
        
        ## Components of the page
        @vfform ?= new GeneForm('Virulence Factor', 'vf')
        @amrform ?= new GeneForm('Antimicrobial Resistance', 'amr')
        @sidebar ?= new Sidebar()

        @selected_vf = m.prop([])
        @selected_amr = m.prop([])
        @selected_genomes = m.prop([])

        return

    controller: (options) =>
        @model()
        @tabCtrl = new mc.Tabs.controller('genes')

        

        
        return


    view: () =>
        ctrl = @tabCtrl 

        tabOptions = {
            flavor: 'bs/nav-tabs'
            tabs: [
                { name: 'genes', label: 'Select Genes'}
                { name: 'genomes', label: 'Select Genomes'}
                { name: 'submit', label: 'Submit Query'}
            ]
        }
 
        renderTabContents = (ctrl) ->
            @activeTab = ctrl.activeTabName
            switch (@activeTab())
                when 'genes' then return renderGeneForm(ctrl)
                when 'genomes' then return renderGenomeSearch(ctrl)
                when 'submit' then return renderSubmit(ctrl)
                else return m('h1', @activeTab())

        renderGeneForm = (ctrl) =>
            return m('div', {class: 'tab-content'}, [
                m('div', {class: 'tab-pane active', id: 'gene-search-querygenes'}, [
                    m('div', {class: 'panel-group genes-search', id: 'accordian'}, [
                        @vfform.view()
                        @amrform.view()
                    ])
                ])
            ])

        renderGenomeSearch = (ctrl) =>
            return m('div', {class: 'tab-content'}, [
                m('div', {class: 'tab-pane active', id: 'gene-search-querygenes'}, [
                    m('div', {class: 'panel-group genes-search', id: 'accordian'}, [
                        "Testing Genome Search"
                    ])
                ])
            ])

        renderSubmit = (ctrl) =>
            ## Setting selected genes/genomes from form
            @selected_vf(@vfform.returnSelected())
            @selected_amr(@amrform.returnSelected())
            @selected_genomes([])

            return m('div', {class: 'tab-content'}, [
                m('div', {class: 'tab-pane active', id: 'gene-search-submit'}, [
                    @vfform.selectedGenesView()
                    @amrform.selectedGenesView()
                    m('div', {class: 'row'}, [
                        m('div', {class: 'col-md-4 col-md-offset-1'}, [
                            m('div', {class: 'panel panel-default'}, [
                                m('div', {id: 'vf-selected-count', class: 'panel-body'}, [
                                    "Selected number of genomes go here"
                                ])
                            ])
                        ])
                    ])
                    m('div', {class: 'row'}, [
                        m('div', {class: 'gene-search-next-wrapper', id: 'query-gene-form'}, [
                            m('button', {class: 'btn btn-success', type: 'submit', value: 'Submit'}, "Submit")
                            m('button', {class: 'btn btn-danger', type: 'reset', value: 'Reset'}, "Reset")
                        ])
                    ])
                ])
            ])

        return [
            super()
            m('div', {id: 'wrapper'}, [
                @sidebar.view()
                m('div', {id: 'page-content-wrapper'}, [
                    m('div', {id: 'page-content -inset'}, [
                        m('div', {class: 'container-fluid'}, [
                            m('div', {class: 'row'}, [
                                m('div', {class: 'col-xs-8'}, [
                                    m('div', {class: 'content-header'}, [
                                        m('h1', [
                                            m('span', {class: "title_part1"}, 'VIRULENCE & AMR')
                                        ])
                                    ])
                                ])
                                m('div', {class: 'col-xs-4'}, [
                                    m('button', {id: "intro-button", class: "btn btn-danger btn-lg", type:"button"}, [
                                        "INTRODUCTION"
                                    ])
                                ])
                            ])
                            m('div', {class: 'intro'}, [
                                m('p', 'BSearch for the presence or absence of virulence factor genes or antimicrobial resistance 
                                        genes in genomes of interest. Detailed information on individual virulence factor or 
                                        antimicrobial resistance genes can be retrieved by clicking on the individual genes.')
                            ])
                            m('.container', [
                                mc.Tabs.view(ctrl, tabOptions)
                                renderTabContents(ctrl)
                            ])
                        ])
                    ])
                ])
            ])
        ]