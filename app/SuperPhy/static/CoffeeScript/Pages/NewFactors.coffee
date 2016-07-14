###
CLASS NewFactors

A page component for the New VF & AMR page

###

# coffeelint: disable=max_line_length

class NewFactors extends Page
    Routes.add("/new", @)
    @controller: (args) ->
        @active = m.prop("genes")
        @tabCtrl = new mc.Tabs.controller('genes')
        @model = GeneSearchModel
        @vfcategories = @model.getCategories('vf')
        @amrcategories = @model.getCategories('amr')
        return @

    @view: (ctrl) ->
        NewrenderTabContents = (ctrl) ->
            activeTab = ctrl.tabCtrl.activeTabName
            switch (activeTab())
                when 'genes' then return NewrenderGeneForm(ctrl)
                when 'genomes' then return NewrenderGenomeSearch(ctrl)
                when 'submit' then return NewrenderSubmit(ctrl)
                else return m('p', "Not a valid tab name.")

        # Below: Defining the 3 functions used above
        NewrenderGeneForm = (ctrl) ->
            m('.tab-content', [
                m('.tab-pane active', {id: 'gene-search-querygenes'}, [
                    m('.panel-group genes-search', {id: 'accordian'}, [
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

        NewrenderGenomeSearch = (ctrl) ->
            m('.tab-content', [
                m('.tab-pane active', {id: 'gene-search-querygenes'}, [
                    m('.panel-group genes-search', {id: 'accordian'}, [
                        "Testing Genome Search (you can change the genomes manually in the GeneSearchModel)"
                    ])
                ])
            ])

        NewrenderSubmit = (ctrl) ->
            m('.tab-content', [
                m('.tab-pane active', {id: 'gene-search-querygenes'}, [
                    m('.panel-group genes-search', id: 'accordian', [
                        m.component(NewSubmitView, {
                            data: ctrl.model
                        })
                    ])
                ])
            ])
        # End function definitions

        super(
            m('.', {id: "wrapper"}, [
                #m.component(Sidebar)
                m('.', {id: 'page-content-wrapper'}, [
                    m('.', {id: 'page-content -inset'}, [
                        m('.container-fluid', [
                            m.component(NewFactorsIntro)
                            m('.container', [
                                mc.Tabs.view(ctrl.tabCtrl,
                                    flavor: 'bs/nav-tabs'
                                    tabs: [
                                        {name: 'genes', label: 'Select Genes'}
                                        {name: 'genomes', label: 'Select Genomes'}
                                        {name: 'submit', label: 'Submit Query'}
                                    ]
                                )
                                NewrenderTabContents(ctrl)
                            ])
                        ])
                    ])
                ])
            ])
        )



###
COMPONENT NewFactorsIntro

Introduction for gene indentification feature

Args:
    none
###
NewFactorsIntro =
    view: (ctrl, args) ->
        m('.intro', [
            m('.row', [
                m('.col-xs-8', [
                    m('.content-header', [
                        m('h1', [
                            m('span.title_part1', 'VIRULENCE & AMR '),
                            m('span.title_part2', 'GENES')
                        ])
                    ])
                ])
                m('.col-xs-4', [
                    m('button.btn btn-danger btn-lg', {id: "intro-button", \
                                                       type:"button"}, [
                        "INTRODUCTION"
                    ])
                ])
            ])
            m('p', "Search for the presence or absence of
                    virulence factor genes or antimicrobial
                    resistance genes in genomes of interest.
                    Detailed information on individual virulence
                    factor or antimicrobial resistance genes
                    can be retrieved by clicking on the
                    individual genes.")
        ])




###
COMPONENT NewSubmitView

View for the third tab on the page 
(This could be moved to another file?)

Args:
    data: gene search model
###
NewSubmitView =
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
                m.component(NewSubmitSelectedView, {
                    selected: args.data.selectedVF
                    title: "Virulence Factor"
                })
                m.component(NewSubmitSelectedView, {
                    selected: args.data.selectedAMR
                    title: "Antimicrobial Resistance"
                })
                m.component(NewSubmitSelectedView, {
                    selected: args.data.selectedGenomes
                    title: "Genome"
                })
                m('.row', [
                    m('.gene-search-next-wrapper', {id: 'query-gene-form'}, [
                        m('button', {class: 'btn btn-success', \
                                     type: 'submit', \
                                     value: 'Submit', \
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

###
COMPONENT NewSubmitSelectedView

Component that confirms number of gene selected

Args:
    selected: list of selected genes or genomes
    title: type of selected genes or genomes
###
NewSubmitSelectedView =
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






