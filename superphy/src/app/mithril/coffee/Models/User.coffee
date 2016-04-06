# coffeelint: disable=max_line_length

Login = ["Username", "Password"]


class FormGroup
    @controller: (args) ->
        args = args or {}
        @placeholder = args.placeholder or ""
        @type = args.type or "text"
        @id = args.id or ""
        @value = m.prop("")
        @help_block = args.help_block or null
        @name = args.name or ""
        @maxlength = args.maxlength or m.prop(20)
        @tabindex = args.tabindex or "1"
        return @
                
    @view: (ctrl) ->
        m('.form-group'
            m('label.col-sm-2 control-label', {for: ctrl.id}, ctrl.placeholder)
                m('.col-sm-8'
                    m('input.form-control', {
                        id:ctrl.id
                        type:ctrl.type
                        tabindex: ctrl.tabindex
                        oninput:m.withAttr('value', ctrl.value)
                        value:ctrl.value()
                        maxlength: ctrl.maxlength
                        size:"15"
                        name:ctrl.name
                        placeholder: ctrl.placeholder
                    })
                    if ctrl.help_block
                        m('span.help-block', "(#{ctrl.help_block})")
                )
        )

class CreateAccount
    @controller: (args) ->
    @view: (ctrl) ->
        m('.container'
            m('.content-header', m('h1'
                m('span.title_part1', "CREATE"),
                m('span.title_part2', "ACCOUNT")
            ))
            m('.panel panel-default', m('panel-body', m('form.form-horizontal',
                {method:"post", action: -> alert("POSTED")}
                m('.form-group', m('.col-sm-offset-2 col-sm-8'
                    m('p', "Already registered or want to make changes to your existing account?",
                        m('a', {tabindex:"9", href:"edit_account"}, "Sign in")
                    )
                    m('p', "To create a new account, fill in all fields below.")
                ))
                m.component(FormGroup, {
                    id:"inputUsername"
                    type:"text"
                    name:"u_username"
                    tabindex:"1"
                    maxlength:"20"
                    placeholder:"Username"
                    help_block:"alphanumeric"
                })
                m.component(FormGroup, {
                    id:"inputPassword"
                    type:"Password"
                    name:"u_password"
                    tabindex:"2"
                    maxlength:"10"
                    placeholder:"Password"
                    help_block:"(6-10 characters)"
                })
                m.component(FormGroup, {
                    id:"inputPassword2"
                    type:"Password"
                    name:"password_confirm"
                    tabindex:"3"
                    maxlength:"10"
                    placeholder:"Re-enter Password"
                })
                m.component(FormGroup, {
                    id:"inputFirstname"
                    type:"text"
                    name:"u_first_name"
                    tabindex:"4"
                    maxlength:"30"
                    placeholder:"First Name"
                })
                m.component(FormGroup, {
                    id:"inputLastname"
                    type:"text"
                    name:"u_last_name"
                    tabindex:"5"
                    maxlength:"30"
                    placeholder:"Last Name"
                })
                m.component(FormGroup, {
                    id:"inputEmail"
                    type:"text"
                    name:"u_email"
                    tabindex:"6"
                    maxlength:"40"
                    placeholder:"Email"
                })
                m('.form-group'
                    m('.col-sm-offset-2 col-sm-8'
                        m('input.btn btn-primary', {
                            type:"submit"
                            tabindex:"7"
                            value:"CreateAccount"
                        })
                    )
                )
            )))
        )