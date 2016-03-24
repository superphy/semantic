class GeneResults extends PageTemplate
    constructor: () ->
        @model()
        @controller()

    model: () =>
        @title = "Virulence Factor and AMR Results"

    controller: () =>

    view: () =>
        return [
            super()
            m('div', {class: 'toc'}, [
                m('div', {class: 'well toc-well'}, [
                    m('div', {id: 'vf_result_legend', class: 'legend', display: 'none'}, [
                        @title
                    ])
                ])
            ])
            m('div', {id: 'results'}, [
                m('div', {id: 'vf_results'},[
                    m('hr')
                    m('h4', "Detected Virulence Factor Alleles")
                ])
            ])
        ]

