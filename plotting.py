def plot_mean(mesh, input_currents, t)
    plt.gca().set_aspect('equal')
    plt.ylim([mesh.boundary_box[1]-5, mesh.boundary_box[3]+20])
    plt.xlim([mesh.boundary_box[0]-5, mesh.boundary_box[2]+5])
    rows = 2
    if opt == 'msfn' and len(functions) == 3:
        rows = 3
    elif len(functions) == 6:
        rows = 4
    #plt.subplot(rows, 2, 1)
    mesh.plot()
    for i, c in enumerate(input_currents):
        mesh.plot_curve(c, color=colors.next(), label='%s, %d'%(functions[i], comb[i]), linewidth=5)
    title = '%s, lambda=%.04f, %s' % (opt, l, str(comb))
    mesh.plot_curve(x, title)
    plt.legend(loc='upper right')

def plot_decomposition(mesh, input_currents, t, q, r, save=True)
    for i, r_i in enumerate(r):
        color = colors.next()
        fig.clf()
        plt.gca().set_aspect('equal')
        plt.ylim([mesh.boundary_box[1]-5, mesh.boundary_box[3]+20])
        plt.xlim([mesh.boundary_box[0]-5, mesh.boundary_box[2]+5])
        mesh.plot()
        mesh.plot_simplices(r_i, color=color)
        mesh.plot_curve(q[i], title=title + ', Q%d&R%d'%(i+1,i+1), color='m', marker='*', linewidth=6, label='Q%d'%(i+1))
        mesh.plot_curve(t, linewidth=4, label="Mean")
        if opt =='msfn' and i== r.shape[0]-1:
            pass
        else:
            mesh.plot_curve(input_currents[i], color='r', ls='--', \
            label='%s, %d'%(functions[i], comb[i]))
        plt.legend(loc='upper right')
        figname = '/home/altaa/fig_dump/%d-%s-%.04f.png'%(figcount,opt,l)
        plt.savefig(figname, dpi=fig.dpi)
        pdf_file.savefig(fig)
        figcount += 1

