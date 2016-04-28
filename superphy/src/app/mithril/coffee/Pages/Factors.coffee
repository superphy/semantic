###
CLASS Factors

Page component for the Factors (VF and AMR) page
###

class Factors extends Page
    Routes.add('/factors', @)
    @controller: (args) ->
        @active = m.prop("genes")
        @tabCtrl = new mc.Tabs.controller('genes')
        @model = GeneSearchModel
        @vfcategories = @model.getCategories('vf')
        @amrcategories= @model.getCategories('amr')
        return @

    @view: (ctrl) ->
        tabOptions = {
            flavor: 'bs/nav-tabs'
            tabs: [
                { name: 'genes', label: 'Select Genes'}
                { name: 'genomes', label: 'Select Genomes'}
                { name: 'submit', label: 'Submit Query'}
            ]
        }

        renderTabContents = (ctrl) ->
            activeTab = ctrl.tabCtrl.activeTabName
            switch (activeTab())
                when 'genes' then return renderGeneForm(ctrl)
                when 'genomes' then return renderGenomeSearch(ctrl)
                when 'submit' then return renderSubmit(ctrl)
                else return m('p', "Not a valid tab name.")

        renderGeneForm = (ctrl) ->
            m('.', {class: 'tab-content'}, [
                m('.', {class: 'tab-pane active', id: 'gene-search-querygenes'}, [
                    m('.', {class: 'panel-group genes-search', id: 'accordian'}, [
                        m.component(GeneSearchPanel, {
                            title: "Virulence Factor"
                            type: "vf"
                            data: ctrl.model.vfList
                            categories: ctrl.vfcategories()
                        })
                        m.component(GeneSearchPanel, {
                            title: "Antimicrobial Resistance"
                            type: "amr"
                            data: ctrl.model.amrList
                            categories: ctrl.vfcategories()
                        })
                    ])
                ])
            ])

        renderGenomeSearch = (ctrl) ->
            m('.', {class: 'tab-content'}, [
                m('.', {class: 'tab-pane active', id: 'gene-search-querygenes'}, [
                    m('.', {class: 'panel-group genes-search', id: 'accordian'}, [
                        "Testing Genome Search (you can change the genomes manually in the GeneSearchModel)"
                    ])
                ])
            ])

        renderSubmit = (ctrl) ->
            m('.', {class: 'tab-content'}, [
                m('.', {class: 'tab-pane active', id: 'gene-search-querygenes'}, [
                    m('.', {class: 'panel-group genes-search', id: 'accordian'}, [
                        m.component(SubmitView, {
                            data: ctrl.model
                        })
                    ])
                ])
            ])


        super(
            m('.', {id: "wrapper"}, [
                #m.component(Sidebar)
                m('.', {id: 'page-content-wrapper'}, [
                    m('.', {id: 'page-content -inset'}, [
                        m('.container-fluid', [
                            m.component(FactorsIntro)
                            m('.container', [
                                mc.Tabs.view(ctrl.tabCtrl, tabOptions)
                                renderTabContents(ctrl)
                            ])
                        ])
                    ])
                ])
            ])
        )

FactorsIntro =
    view: (ctrl, args) ->
        m(".intro", [
            m('.', {class: 'row'}, [
                m('.col-xs-8', [
                    m('.content-header', [
                        m('h1', [
                            m('span', {class: "title_part1"}, 'VIRULENCE & AMR '),
                            m('span', {class: "title_part2"}, 'GENES')
                        ])
                    ])
                ])
                m('.', {class: 'col-xs-4'}, [
                    m('button', {id: "intro-button", \
                                 class: "btn btn-danger btn-lg", \
                                 type:"button"}, [
                        "INTRODUCTION"
                    ])
                ])
            ])
            m("p", "Search for the presence or absence of
                    virulence factor genes or antimicrobial
                    resistance genes in genomes of interest.
                    Detailed information on individual virulence
                    factor or antimicrobial resistance genes
                    can be retrieved by clicking on the
                    individual genes.")
        ])

SubmitView =
    controller: (args) ->
        @submit = () ->
            args.data.submit()
        @reset = () ->
            args.data.reset()
        return @

    view: (ctrl, args) ->
        args.data.getSelectedVF()
        args.data.getSelectedAMR()
        m('.tab-content', [
            m('.tab-pane active', {id: 'gene-search-submit'}, [
                m.component(SubmitSelectedView, {
                    selected: args.data.selectedVF
                    title: "Virulence Factor"
                })
                m.component(SubmitSelectedView, {
                    selected: args.data.selectedAMR
                    title: "Antimicrobial Resistance"
                })
                m.component(SubmitSelectedView, {
                    selected: args.data.selectedGenomes
                    title: "Genome"
                })
                m('.row', [
                    m('.gene-search-next-wrapper', {id: 'query-gene-form'}, [
                        m('button', {class: 'btn btn-success', \
                                     type: 'submit', \
                                     value: 'Submit',
                                     onclick: ctrl.submit},
                           "Submit")
                        m('button', {class: 'btn btn-danger', \
                                     type: 'reset', \
                                     value: 'Reset', \
                                     onclick: ctrl.reset}, 
                           "Reset")
                    ])
                ])
            ])
        ])

SubmitSelectedView =
    view: (ctrl, args) ->
        m('.row', [
            m('.col-md-4 col-md-offset-1', [
                m('.panel panel-default', [
                    m('.panel-body', {id: 'vf-selected-count'}, [
                        args.selected.length
                        " "
                        args.title.toLowerCase()
                        if args.selected.length isnt 1 then "s"
                        " selected"
                    ])
                ])
            ])
        ])






