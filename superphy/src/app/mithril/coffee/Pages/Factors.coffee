class Factors extends PageTemplate
    #a list to be passed to the link generator function

    model: () =>
        @value = m.prop('')
        @activeTab = m.prop('genes')
        
        
        @vfform ?= new GeneForm('Virulence Factor', 'vf')
        @amrform ?= new GeneForm('Antimicrobial Resistance', 'amr')
        @sidebar ?= new Sidebar()

        return

    controller: (options) =>
        @model()
        @tabCtrl = new mc.Tabs.controller('genes')
        
        return


    view: () =>
        ctrl = @tabCtrl # Need to fix this so that it goes back within the controller or model @___@
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
            return m('div', {class: 'tab-content'}, [
                m('div', {class: 'tab-pane active', id: 'gene-search-submit'}, [
                    m('div', {class: 'row'}, [
                        m('div', {class: 'col-md-4 col-md-offset-1'}, [
                            m('div', {class: 'panel panel-default'}, [
                                m('div', {id: 'vf-selected-count', class: 'panel-body'}, [
                                    "Selected number of genomes go here"
                                ])
                            ])
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
                ])
            ])
        ]