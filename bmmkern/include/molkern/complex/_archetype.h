#ifndef _ARCHETYPE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H
#define _ARCHETYPE__F9ED1116_EDB9_5e17_25FF_F745B15D0100__H

#include "molkern/__moldefs.h"

#include "molkern/forcefield/_forcefield.h"
#include "molkern/forcefield/_forcefield_amber.h"
#include "molkern/forcefield/_residome.h"
#include "molkern/forcefield/_residome_amber.h"
#include "molkern/forcefield/_bond.h"
#include "molkern/forcefield/_angle.h"
#include "molkern/forcefield/_torsion.h"
#include "molkern/forcefield/_rotamer.h"

#include "molkern/complex/_region.h"
#include "molkern/complex/_geom_tool.h"
#include "molkern/complex/_protonization.h"

namespace molkern
{
  using namespace prgkern;

  /**
  *  Это образец, (шаблон, архетип, ..) с которого снимаются копии молекул.
  *  Он позволяет хранить данные единые для всех молекул заданного типа, и
  *  вычислять единые величины (энергии связей, углов, и т.д.). Результаты
  *  вычислений, естественно, для разных молекул будут разными, поскольку
  *  каждая конкретная молекула заданного типа имеет разное число свободных
  *  атомов и различное положение атомов в пространстве.
  */

  /**
  *  Шаблон для архетипов молекул, считываемых в любых силовых полях.
  */
  template
  <
    int FORCEFIELD_TYPE = FORCEFIELD_AMBER_,
    int RESIDOME_TYPE = RESIDOME_AMBER_
  >
  class Archetype_;

  /**
  *  Шаблон архетипа для силового поля AMBER + GAFF.
  *  Данный архетип позволяет работать только со следующими компонентами
  *  молекулы: связи, валентные углы, торсионные углы, 1-4 взаимодействия,
  *  связанные цепи атомов (актуально для растворов), ротамеры.
  */
  template <>
  class Archetype_<FORCEFIELD_AMBER_, RESIDOME_AMBER_>
  {
    typedef Residome_  <RESIDOME_AMBER_>    _Residome;
    typedef Forcefield_<FORCEFIELD_AMBER_>  _Forcefield;
    typedef _Residome::residue_type         _Residue;

    typedef Atomdata_                       _Atomdata;
    typedef Bond_<FORCEFIELD_AMBER_>        _Bond;
    typedef Angle_<FORCEFIELD_AMBER_>       _Angle;
    typedef Torsion_<FORCEFIELD_AMBER_>     _Torsion;
    typedef Interaction                     _Interaction;

    typedef index_<2>                       _Pair14;
    typedef Rotamer_                        _Rotamer;
    typedef RotamerConnect_                 _Edge;
    typedef Chain_                          _Chain;
    typedef Box_<3, real_t>                 _Box;

  public:

    enum { dimension = 3 };

    typedef _Forcefield   forcefield_type;
    typedef _Residome     residome_type;
    typedef _Atomdata     atomdata_type;
    typedef _Atomdata     atom_type;
    typedef _Bond         bond_type;
    typedef _Angle        angle_type;
    typedef _Torsion      torsion_type;
    typedef _Rotamer      rotamer_type;
    typedef _Chain        chain_type;
    typedef _Edge         edge_type;

    /**
     * @param forcefield - база данных по силовому полю
     * @param residome - база данных по топологии
     * @param freedom_type - предельный тип свободы, допустимый для всех молекул архетипа
     */
    Archetype_(const _Forcefield *forcefield, const _Residome *residome, unsigned freedom_type)
    : forcefield_(forcefield), residome_(residome), freedom_type_(freedom_type),
      filename_(""), molname_(""), format_(FORMAT_UNKNOWN_), is_solution_(false)
    {}

    /**
    * @brief загружает и строит прототип молекулы из файла
    * @note формат файла определяется расширением имени файла
    * @param filename полное (с путем) имя файла
    * @param altpos выбор альтернативной позиции (для *.pdb)
    *   Заметим, что не всегда лучшей позицией является 'A'.
    */
    bool load(const std::string &filename, char altpos='A');
    bool load(_I2T<FORMAT_PDB_>, std::ifstream &file, char altpos='A');
    bool load(_I2T<FORMAT_HIN_>, std::ifstream &file, char altpos='A');
    bool load(_I2T<FORMAT_MOL2_>,std::ifstream &file, char altpos='A');
    bool load(_I2T<FORMAT_BMM_>, std::ifstream &file, char altpos='A');
    bool load(_I2T<WATER_>, const std::string &solution_name);

    /**
    * @brief построение недостающих взаимосвязей для того или иного формата после загрузки
    * @param pH кислотность воды, для которой строятся компоненты (действительно только для PDB)
    */
    void build(_I2T<MOLECULE_>, int pH=NORMAL_PH_WATER);
    bool build(_I2T<WATER_>,    int pH=NORMAL_PH_WATER);

    /**
    * @brief сохраняем молекулу в файл, причем формат определяется расширением
    * @param filename полное имя файла (с путем)
    * @param atom стартовый атом в базе данных атомов
    * @param prn_hydrogens печатать (или нет) атомы водорода
    */
    template <typename _Atom> unsigned save(const std::string &filename,
      const _Atom *atoms, bool prn_hydrogens=YES_HYDROGEN) const;
    template <typename _Atom> unsigned save(_I2T<FORMAT_PDB_ >, std::ofstream &file,
      const _Atom *atoms, bool prn_hydrogens=YES_HYDROGEN) const;
    template <typename _Atom> unsigned save(_I2T<FORMAT_HIN_ >, std::ofstream &file,
      const _Atom *atoms, bool) const;
    template <typename _Atom> unsigned save(_I2T<FORMAT_MOL2_>, std::ofstream &file,
      const _Atom *atoms,  bool) const;
    template <typename _Atom> unsigned save(_I2T<FORMAT_BMM_ >, std::ofstream &file,
      const _Atom *atoms, bool) const;

    /**   ФУНКЦИИ ДЛЯ НАХОЖДЕНИЯ ОБЪЕКТОВ, СВЯЗАННЫХ С ЗАДАННЫМИ АТОМАМИ
    *  Данные функции используют, например, для получения подмножетва связей,
    *  углов и т.д., которые связаны со списком свободных атомов. Имея такое
    *  подмножество, при расчетах энергий можно не вычислять вклады, которые
    *  связаны с фиксированными атомами, и которые дают всего лишь константу
    *  в полную энергию. В силу самоочевидности этих функций подробное
    *  документирование не используется.
    *
    * @param start,end итераторы последовательности индексов
    * @param to итератор начала вывода данных
    * @return количество выведенных данных
    *
    *  Итераторы определяют индексы, по которым из соответствующих масивов
    *  выдергиваются соответствующие элементы. Так можно проходить, как и
    *  выделенные элементы (связанные со свободными атомами), так и все элементы.
    *  Для последнего случая нужно задать итераторы со значениями 0 и size().

    * @note (!)Пользователь должен сам заботиться о том, чтобы в буфере вывода
    *  было достаточно места. Если это не так, то резуьтат вывода неопределен.
    *  В выводе идут беззнаковые числа, соответствующие номерам связей,
    *  которые включают хотя бы один атом из заданной последовательности
    */

    template <typename __Iterator, typename _Iterator>
    unsigned free(_I2T<BOND_>, __Iterator &to, _Iterator it, _Iterator ite) const;

    template <typename __Iterator, typename _Iterator>
    unsigned free(_I2T<ANGLE_>, __Iterator &to, _Iterator it, _Iterator ite) const;

    template <typename __Iterator, typename _Iterator>
    unsigned free(_I2T<TORSION_>, __Iterator &to, _Iterator it, _Iterator ite) const;

    template <typename __Iterator, typename _Iterator>
    unsigned free(_I2T<PAIR14_>, __Iterator &to, _Iterator it, _Iterator ite) const;

    /**       ФУНКЦИИ ДЛЯ РАСЧЕТА ПОТЕНЦИАЛЬНОЙ ЭНЕРГИИ И ЕЕ ПРОИЗВОДНЫХ
    *  Данные функции используют для получения разнообразных энергетических
    *  вкладов в потенциальную энергию молекулы и для расчета ее производных.
    *  В силу самоочевидности этих функций подробное документирование
    *  не используется.
    *
    * @param atoms стартовый элемент внешнего массива атомов
    * @param start,end итераторы последовательности индексов
    *
    *  Итераторы определяют индексы, по которым из соответствующих масивов
    *  выдергиваются соответствующие элементы. Так можно проходить, как и
    *  выделенные элементы (связанные со свободными атомами), так и все элементы.
    *  Для последнего случая нужно задать итераторы со значениями 0 и array.size.
    *  Заметим, чтобы пройти по элементам, связанными со свободными атомами,
    *  необходимо сперва вызвать функции build(_I2T<FREE_..>, ..), которые
    *  определяют номера этих элементов.
    */

    template <typename _Atom, typename _Iterator>
    _E(real_t) U(_I2T<BOND_>,  const _Atom *atoms, _Iterator start, _Iterator end,
      bool make_print=YES_PRINT) const;

    template <typename _Atom, typename _Iterator>
    _E(real_t) U(_I2T<ANGLE_>, const _Atom *atoms, _Iterator start, _Iterator end,
      bool make_print=YES_PRINT) const;

    template <typename _Atom, typename _Iterator>
    _E(real_t) U(_I2T<TORSION_>, const _Atom *atoms, _Iterator start, _Iterator end,
      bool make_print=YES_PRINT) const;

    template <typename _Atom, typename _Iterator>
    _E(real_t) U(_I2T<PAIR14_>,const _Atom *atoms, _Iterator start, _Iterator end,
      bool make_print=YES_PRINT) const;

    template <typename _Atom, typename _Iterator>
    _E(real_t) dU__dX(_I2T<BOND_>,   _Atom *atoms, _Iterator start, _Iterator end) const;

    template <typename _Atom, typename _Iterator>
    _E(real_t) dU__dX(_I2T<ANGLE_>,  _Atom *atoms, _Iterator start, _Iterator end) const;

    template <typename _Atom, typename _Iterator>
    _E(real_t) dU__dX(_I2T<TORSION_>, _Atom *atoms, _Iterator start, _Iterator end) const;

    template <typename _Atom, typename _Iterator>
    _E(real_t) dU__dX(_I2T<PAIR14_>, _Atom *atoms, _Iterator start, _Iterator end) const;

    /**           ФУНКЦИИ ДЛЯ ПОЛУЧЕНИЯ ИНФОРМАЦИИ О МОЛЕКУЛЕ
    *  Данные функции используют для получения разнообразной информации
    *  о молекуле. В силу самоочевидности этих функций подробное документирование
    *  не используется.
    */
    DEFINE_VECTOR_ACCESS_FUNCTION(ATOMDATA_, _Atomdata, atomdata_);
    DEFINE_VECTOR_ACCESS_FUNCTION(ATOM_,     _Atomdata, atomdata_);
    DEFINE_VECTOR_ACCESS_FUNCTION(BOND_,     _Bond,     bonds_   );
    DEFINE_VECTOR_ACCESS_FUNCTION(ANGLE_,    _Angle,    angles_  );
    DEFINE_VECTOR_ACCESS_FUNCTION(TORSION_,  _Torsion,  torsions_);
    DEFINE_VECTOR_ACCESS_FUNCTION(PAIR14_,   _Pair14,   pair14s_ );
    DEFINE_VECTOR_ACCESS_FUNCTION(ROTAMER_,  _Rotamer,  rotamers_);
    DEFINE_VECTOR_ACCESS_FUNCTION(CHAIN_,    _Chain,    chains_  );
    DEFINE_VECTOR_ACCESS_FUNCTION(EDGE_,     _Edge,     edges_   );
    DEFINE_VECTOR_ACCESS_FUNCTION(ROOT_ROTAMER_, unsigned, roots_);

    const std::string name(_I2T<FILE_>    ) const { return filename_; }
    const std::string name(_I2T<MOLECULE_>) const { return molname_; }
    const std::string name() const
    {
      std::string msg = _S("\"") + filename_ + _S("\"<") + molname_ + _S(">");
      return msg;
    }

    /// максимальный тип свободы архетипа
    unsigned get(_I2T<FREEDOM_TYPE_>) const { return freedom_type_; }

    /// полное число степеней свободы для разных типов движения
    unsigned count(_I2T<FREEDOM_>, unsigned freedom_type) const
    {
      unsigned cnt = 0;

      if (freedom_type & YES_ATOM_)
      {
        if (freedom_type & YES_ROTAMER_)
        {
          cnt += freedom_count_(_I2T<YES_ATOM_ | YES_ROTAMER_>());
        }
        else
        {
          cnt += freedom_count_(_I2T<YES_ATOM_ |  NO_ROTAMER_>());
          if (freedom_type & NO_CM_)
          {
            if (freedom_type & YES_UNION_) cnt -= freedom_count_(_I2T<YES_CM_ | YES_UNION_>());
            else cnt -= freedom_count_(_I2T<YES_CM_ |  NO_UNION_>());
          }
        return cnt;
        }
      }

      if (freedom_type & YES_ROTAMER_) cnt += freedom_count_(_I2T<YES_ROTAMER_>());

      if (freedom_type & YES_CM_)
      {
        if (freedom_type & YES_UNION_) cnt += freedom_count_(_I2T<YES_CM_ | YES_UNION_>());
        else cnt += freedom_count_(_I2T<YES_CM_ |  NO_UNION_>());
      }

      return cnt;
    }

    /**
    *  Рассчитывает приблизительный ящик, к котором могут разместиться все атомы
    *  молекулы. Оценка дается по самым отклоняющимся атомам.
    * @return приблизительный ящик, в котором могут разместится все атомы молекулы
    */
    _Box get(_I2T<BOX_>) const
    {
      unsigned atom_count = count(_I2T<ATOM_>());
      return calculate(BOX, &atomdata_[0], range_iterator(0), range_iterator(atom_count));
    }

    /// проверка, являтся ли молекула раствором
    bool is_solution() const { return is_solution_; }

    /**
    *  Конвертирует итератор по индексам в итератор по Atomdata молекулы.
    *  Позволяет пользователю, который создал группы атомов (через совокупность
    *  индексов) снаружи класса молекулы, получать соответствующие атомы группы.

    *  (!) Разрешаем только константный итератор, чтобы запретить случайное
    *  изменение содержимого atomdata.
    * @param it итератор (начала или конца) последовательности индексов
    * @return итератор по атомам молекулы
    */
    template <typename _Iterator>
    const_array_iterator<_Atomdata, _Iterator> make_iterator(_Iterator it) const
    {
      return const_array_iterator<_Atomdata, _Iterator>(&atomdata_[0], it);
    }

    /**
     *   запись обобщенной позиции по заданному адресу
     * @param x - адрес записи
     * @param X - координаты атомов
     * @param angle - угол поворота относительно atomdata
     * @param ranx - смещение относительно atomdata
     * @return число записанных позиций
     */
    template <typename _Atom> unsigned
    write_position(_I2T<MOLECULE_>, real_t *x, unsigned freedom_type, const _Atom *atoms,
      const vector_t &angle, const vector_t &ranx) const
    {
      unsigned cnt = 0;
      if (freedom_type & YES_ATOM_) cnt += write_position_(_I2T<YES_ATOM_>(), x + cnt, atoms);
      if (freedom_type & YES_ROTAMER_) cnt += write_position_(_I2T<YES_ROTAMER_>(), x + cnt, atoms);
      if (freedom_type & YES_CM_)
        cnt += write_position_(_I2T<YES_CM_ | YES_UNION_>(), x + cnt, atoms, angle, ranx);
        // начальная позиция формируется одинаково для UNION и без UNION
      return cnt;
    }

    template <typename _Atom> unsigned
    write_position(_I2T<WATER_>, real_t *x, unsigned freedom_type, const _Atom *atoms) const
    {
      if (freedom_type & YES_ATOM_) return write_position_(_I2T<YES_ATOM_>(), x, atoms);

      unsigned cnt = 0;
      if (freedom_type & YES_CM_)
      {
        cnt += 6 * count(CHAIN);
        for (unsigned i=0; i<cnt; i++) x[i] = 0.f;
      }
      return cnt;
    }

    template <typename _Atom> unsigned
    read(_I2T<POSITION_>, const real_t *x, unsigned freedom_type, _Atom *atoms) const
    {
      unsigned cnt = 0; bool use_atomdata = true;

      if (freedom_type & YES_ATOM_)
        cnt += read_position_(_I2T<YES_ATOM_>(), x + cnt, atoms, use_atomdata);

      if (freedom_type & YES_ROTAMER_)
        cnt += read_position_(_I2T<YES_ROTAMER_>(), x, atoms, use_atomdata);
      if (freedom_type & YES_CM_)
      {
        if (freedom_type & YES_UNION_)
          cnt += read_position_(_I2T<YES_CM_ | YES_UNION_>(), x + cnt, atoms, use_atomdata);
        else cnt += read_position_(_I2T<YES_CM_ | NO_UNION_>(), x + cnt, atoms, use_atomdata);
      }
      return cnt;
    }

    template <typename _Atom> unsigned
    write(_I2T<GRADIENT_>, real_t *g, unsigned freedom_type, const _Atom *atoms) const
    {
      unsigned cnt = 0;
      if (freedom_type & YES_ATOM_) cnt += write_gradient_(_I2T<YES_ATOM_>(), g + cnt, atoms);
      if (freedom_type & YES_ROTAMER_) cnt += write_gradient_(_I2T<YES_ROTAMER_>(), g, atoms);
      if (freedom_type & YES_CM_)
      {
        if (freedom_type & YES_UNION_)
          cnt += write_gradient_(_I2T<YES_CM_ | YES_UNION_>(), g + cnt, atoms);
        else cnt += write_gradient_(_I2T<YES_CM_ | NO_UNION_>(), g + cnt, atoms);
      }
      return cnt;
    }

    template <typename _Atom> unsigned
    read_position_(_I2T<YES_ATOM_>, const real_t *x, _Atom *atoms, bool &use_atomdata) const
    {
      unsigned atom_count = count(ATOM);
      for (unsigned i=0, k=0; i<atom_count; i++, k+=3)
      {
        _Atom &atom = atoms[i];
        atom.X[0] = x[k    ];
        atom.X[1] = x[k + 1];
        atom.X[2] = x[k + 2];
      }
      use_atomdata = false;

      return 3 * atom_count;
    }

    template <typename _Atom> unsigned
    read_position_(_I2T<YES_ROTAMER_>, const real_t *x, _Atom *atoms, bool &use_atomdata) const
    {
      // Порядок обобщенных координат следующий: (1) углы вращения всех ротамеров;
      // (2) углы вращения сформированных цепей; (3) декартовые смещения сформированных цепей.
      // Порядок обусловлен тем, что не меняется центр тяжести молекулы после установки углов
      // ротамеров. Использование вместо Эйлеровых углов вращения вокруг осей x, y, z, позволяет
      // достаточно просто подсчитывать момент вращения всей молекулы.

      //-----------------------------------------------------------------------------
      //                       перекопируем корни
      //-----------------------------------------------------------------------------
      const unsigned *roots = get(ROOT_ROTAMER);
      unsigned root_count = count(ROOT_ROTAMER);
      const _Atomdata *atomdata = get(ATOMDATA);

      for (unsigned iroot=0; iroot<root_count; iroot++)
      {
        // копируем атомы тела ротамера
        const _Rotamer *rotamer = get(ROTAMER, roots[iroot]);
        const unsigned *ndx = rotamer->get(ATOM);
        unsigned count = rotamer->count(ATOM);
        for (unsigned i=0; i<count; i++)
          atoms[ndx[i]].X = atomdata[ndx[i]].X;

        // копируем атомы стиков, так как они нужны для привязки
        unsigned stick_count = rotamer->count(STICK);
        for (unsigned i=0; i<stick_count; i++)
        {
          unsigned ne = rotamer->get(ROTAMER_EXT_ATOM, i);
          atoms[ne].X = atomdata[ne].X;
        }
      }
      use_atomdata = false;

      //-----------------------------------------------------------------------------
      //                       склеим ротамеры друг с другом
      //-----------------------------------------------------------------------------
      unsigned edge_count = count(EDGE);
      if (edge_count == 0) return 0;

      const _Edge *edges = get(EDGE);
      const _Rotamer *rotamers = get(ROTAMER);

      for (unsigned iedge=0; iedge<edge_count; iedge++)
      {
        const _Edge &edge = edges[iedge];
        unsigned r1 = edge.rotamer_from;
        unsigned r2 = edge.rotamer_to;
        rotamer_paste(atoms, atomdata, rotamers[r1], rotamers[r2],
          edge.stick_from, edge.stick_to, x[iedge]);
      }

      return edge_count;
    }

    template <typename _Atom> unsigned
    read_position_(_I2T<YES_CM_ | NO_UNION_>, const real_t *x, _Atom *atoms, bool &use_atomdata) const
    {
      // считываем данные из хранилища
      if (use_atomdata)
      {
        unsigned atom_count = count(ATOM);
        for (unsigned i=0; i<atom_count; i++) atoms[i].X = atomdata_[i].X;
        use_atomdata = false;
      }

      const _Chain *chains = get(CHAIN);
      unsigned chain_count = count(CHAIN);
      for (unsigned ichain=0; ichain<chain_count; ichain++, x+=6)
      {
        vector_t alpha = vector_t(x[0], x[1], x[2]);
        vector_t X     = vector_t(x[3], x[4], x[5]);
        const _Chain &chain = chains[ichain];
        const unsigned *ndx = chain.get(ATOM);
        unsigned count = chain.count(ATOM);

        vector_t cm = calculate(MASS_CENTER, atoms, ndx, ndx + count);
        rotate(XYZ_ROTATOR, atoms, alpha, cm, ndx, ndx + count);
        move(atoms, X, ndx, ndx + count);
      }
      return 6 * chain_count;
    }

    template <typename _Atom> unsigned
    read_position_(_I2T<YES_CM_ | YES_UNION_>, const real_t *x, _Atom *atoms, bool &use_atomdata) const
    {
      // считываем данные из хранилища
      if (use_atomdata)
      {
        unsigned atom_count = count(ATOM);
        for (unsigned i=0; i<atom_count; i++) atoms[i].X = atomdata_[i].X;
        use_atomdata = false;
      }

      unsigned atom_count = count(ATOM);
      vector_t alpha = vector_t(x[0], x[1], x[2]);
      vector_t X     = vector_t(x[3], x[4], x[5]);

      vector_t cm = calculate(MASS_CENTER, atoms, range_iterator(0), range_iterator(atom_count));
      rotate(XYZ_ROTATOR, atoms, alpha, cm, range_iterator(0), range_iterator(atom_count));
      move(atoms, X, range_iterator(0), range_iterator(atom_count));

      return 6;
    }

    template <typename _Atom> unsigned
    write_gradient_(_I2T<YES_ATOM_>, real_t *g, const _Atom *atoms) const
    {
      unsigned atom_count = count(ATOM);
      for (unsigned i=0, k=0; i<atom_count; i++, k+=3)
      {
        const _Atom &atom = atoms[i];
        g[k    ] = -atom.F[0];
        g[k + 1] = -atom.F[1];
        g[k + 2] = -atom.F[2];
      }
      return 3 * atom_count;
    }

    template <typename _Atom> unsigned
    write_gradient_(_I2T<YES_ROTAMER_>, real_t *g, const _Atom *atoms) const
    {
      //-------------------------------------------------------------------------------------
      //             обработаем вращения вокруг ротамеров
      //-------------------------------------------------------------------------------------
      // код полагается на то, что массив ребер является упорядоченным и сгруппированным по
      // несвязанным группам ротамеров

      typedef boost::tuple<vector_t, vector_t, unsigned, unsigned>  _Tuple;
        // осей вращения, опорная точка вращения, ссылки на момент и ротамер

      unsigned cnt = count(ROTAMER)- count(ROOT_ROTAMER);
      for (unsigned i=0; i<cnt; i++) g[i] = 0.;
        // очистка от старых значений, так как используем накопление в этот массив

      unsigned edge_count = count(EDGE);
      _Tuple tuples[edge_count]; // массив осей вращения, ссылок на моменты и ротамеры

      unsigned cur_tuple = 0;
      for (unsigned iedge=0; iedge<edge_count; iedge++)
      {
        // определяем характеристики текущего ребра
        const _Edge *edge = get(EDGE, iedge);
        unsigned prev_rotamer = edge->rotamer_from;
        unsigned r2 = edge->rotamer_to;
        unsigned s2 = edge->stick_to;

        // любое ребро вводит новую ось вращения, получим ее характеристики
        const _Rotamer *rotamer = get(ROTAMER, r2);
        unsigned t1 = rotamer->get(ROTAMER_EXT_ATOM, s2);
        unsigned t2 = rotamer->get(ROTAMER_INT_ATOM, s2);
        const _Atom &atom1 = atoms[t1];
        const _Atom &atom2 = atoms[t2];
        vector_t W = atom2.X - atom1.X;
        W.normalize(); // так как функции работают с нормализованными векторами для скорости

        // определяем позицию приклеивания нового тапла
        while (cur_tuple > 0)
        {
          _Tuple &tuple = tuples[cur_tuple - 1];
          unsigned r2 = tuple.get<3>();
          if (r2 == prev_rotamer) break;
          cur_tuple--;
        }

        _Tuple &tuple = tuples[cur_tuple++];
        tuple.get<0>() = W;
        tuple.get<1>() = atom2.X;
        tuple.get<2>() = iedge;
        tuple.get<3>() = r2;

        // считаем момент ротамера относительно текущей и предыдущих осей
        for (unsigned im = 0; im<cur_tuple; im++)
        {
          _Tuple &tuple__ = tuples[im];
          vector_t W__ = tuple__.get<0>();
          vector_t O__ = tuple__.get<1>();
          unsigned m__ = tuple__.get<2>();
          g[m__] -= rotamer_moment(atoms, rotamer, W__, O__);
        }
      }
      return edge_count;
    }

    template <typename _Atom> unsigned
    write_gradient_(_I2T<YES_CM_ | NO_UNION_>, real_t *g, const _Atom *atoms) const
    {
      const _Chain *chains = get(CHAIN);
      unsigned chain_count = count(CHAIN);
      for (unsigned ichain=0; ichain<chain_count; ichain++, g+=6)
      {
        real_t mx = 0., my = 0., mz = 0.;
        vector_t F(0., 0., 0.);

        const _Chain &chain = chains[ichain];
        const unsigned *ndx = chain.get(ATOM);
        unsigned count = chain.count(ATOM);
        vector_t cm = calculate(MASS_CENTER, atoms, ndx, ndx + count);
        for (unsigned i=0; i<count; i++)
        {
          const _Atom &atom = atoms[ndx[i]];
          vector_t F__ = atom.F;
          mx += atom_moment(vector_t(1., 0., 0.), cm, F__, atom.X);
          my += atom_moment(vector_t(0., 1., 0.), cm, F__, atom.X);
          mz += atom_moment(vector_t(0., 0., 1.), cm, F__, atom.X);
          F += F__;
        }
        g[0] = -mx;
        g[1] = -my;
        g[2] = -mz;
        g[3] = -F[0];
        g[4] = -F[1];
        g[5] = -F[2];
      }
      return 6 * chain_count;
    }

    template <typename _Atom> unsigned
    write_gradient_(_I2T<YES_CM_ | YES_UNION_>, real_t *g, const _Atom *atoms) const
    {
      unsigned atom_count = count(ATOM);
      real_t mx = 0., my = 0., mz = 0.;
      vector_t F(0., 0., 0.);

      vector_t cm = calculate(MASS_CENTER, atoms, range_iterator(0), range_iterator(atom_count));
      for (unsigned i=0; i<atom_count; i++)
      {
        const _Atom &atom = atoms[i];
        vector_t F__ = atom.F;
        mx += atom_moment(vector_t(1., 0., 0.), cm, F__, atom.X);
        my += atom_moment(vector_t(0., 1., 0.), cm, F__, atom.X);
        mz += atom_moment(vector_t(0., 0., 1.), cm, F__, atom.X);
        F += F__;
      }
      g[0] = -mx;
      g[1] = -my;
      g[2] = -mz;
      g[3] = -F[0];
      g[4] = -F[1];
      g[5] = -F[2];
      return 6;
    }

  protected:

    /**
    * @brief найти атом с данным внешним идентификаторов, если он есть
    * @param sid целевой идентификатор
    * @param pos приблизительная позиция старта поиска
    * @return внутренний номер атома (nill, если не найден)
    */
    int find(_I2T<EXT_ID_>, int sid, int pos=0) const;

    /**
    * @brief конвертировать внутренний идентификатор во внешний
    * @param ndx номер атома
    * @return внешний идентификатор атома (или 0 при его отсутствии)
    */
    int get(_I2T<EXT_ID_>, int internal_ndx) const
    {
      assert(_LT(internal_ndx, (int)atomdata_.size()));
      return atomdata_[internal_ndx].sid;
    }

    /**             ФУНКЦИИ ВОССТАНОВЛЕНИЯ ТОПОЛОГИИ
    *  Данные функции используют допущения, что в pdb файлах все атомы
    *  аминокислотных остатков, которые формируют цепь (N-CA-C-O), определены.
    *  В цепи атомы "держат друг дружку", потому нет неопределенности их
    *  позиций в эксперименте. Исключением являются боковые цепи, которые
    *  могут иметь различное положение в разных белках. В этом случае в pdb
    *  файле возникают разрывы в нумерации аминокислотных остатков.
    *  Функции обрабатывает разрывы путем добавления дополнительных
    *  атомов водорода, поскольку связываение цепей через разрыв приводит
    *  к "дикому" искажению реальной геометрии белка в местах разрывов.
    *  При восстановлении появляются новые атомы для которых отсутствуют
    *  номера для идентификации (sid==0). Предполагается, что нет необходимости
    *  вставлять эти номера. Подобная вставка изменяет идентификаторы атомов,
    *  что затрудняет для пользователя поиск нужных атомов в исходном файле.
    */
    bool build_(_I2T<FORMAT_PDB_  >, int pH=NORMAL_PH_WATER);
    bool build_(_I2T<FORMAT_HIN_  >, int pH=NORMAL_PH_WATER);
    bool build_(_I2T<FORMAT_MOL2_ >, int pH=NORMAL_PH_WATER);
    bool build_(_I2T<FORMAT_BMM_  >, int pH=NORMAL_PH_WATER);
    void build_(_I2T<RESIDUE_NAME_>, int pH=NORMAL_PH_WATER);
    void build_(_I2T<HYDROGEN_       >);
    void build_(_I2T<RESIDUE_CONTACT_>);
    void build_(_I2T<BOND_   >);
    void build_(_I2T<ANGLE_  >);
    void build_(_I2T<TORSION_>);
    void build_(_I2T<PAIR14_ >);
    void build_(_I2T<CONNECT_>);
    void build_(_I2T<ROTAMER_>);
    void build_(_I2T<CHAIN_  >);

    unsigned freedom_count_(_I2T<YES_ROTAMER_>) const { return count(ROTAMER) - count(ROOT_ROTAMER); }
    unsigned freedom_count_(_I2T<YES_CM_ |  NO_UNION_>) const { return 6 * count(CHAIN); }
    unsigned freedom_count_(_I2T<YES_CM_ | YES_UNION_>) const { return 6; }
    unsigned freedom_count_(_I2T<YES_ATOM_ |  NO_ROTAMER_>) const { return 3 * atomdata_.size(); }
    unsigned freedom_count_(_I2T<YES_ATOM_ | YES_ROTAMER_>) const
    {
      unsigned root_count = count(ROOT_ROTAMER);
      const unsigned *ndx = get(ROOT_ROTAMER);
      unsigned atom_count = 0;
      for (unsigned i=0; i<root_count; i++)
      {
        const _Rotamer &rotamer = rotamers_[ndx[i]];
        atom_count += rotamer.count(ATOM);
      }
      return 3 * atom_count;
    }

    template <typename _Atom> unsigned
    write_position_(_I2T<YES_ATOM_>, real_t *x, const _Atom *atoms) const
    {
      unsigned atom_count = count(ATOM);
      for (unsigned i=0, k=0; i<atom_count; i++, k+=3)
      {
        const _Atom &atom = atoms[i];
        x[k    ] = atom.X[0];
        x[k + 1] = atom.X[1];
        x[k + 2] = atom.X[2];
      }
      return 3 * atom_count;
    }

    template <typename _Atom> unsigned
    write_position_(_I2T<YES_ROTAMER_>, real_t *x, const _Atom *atoms) const
    {
      unsigned edge_count = count(EDGE);
      const _Edge *edges_ = get(EDGE);
      const _Rotamer *rotamers_ = get(ROTAMER);
      for (unsigned i=0; i<edge_count; i++)
      {
        const _Edge &edge = edges_[i];
        unsigned r1 = edge.rotamer_from;
        unsigned r2 = edge.rotamer_to;
        unsigned s1 = edge.stick_from;
        unsigned s2 = edge.stick_to;

        const _Rotamer &rotamer1 = rotamers_[r1];
        const _Rotamer &rotamer2 = rotamers_[r2];
        unsigned na = rotamer1.get(ROTAMER_ADD_ATOM, s1);
        unsigned nb = rotamer1.get(ROTAMER_INT_ATOM, s1);
        unsigned nc = rotamer2.get(ROTAMER_INT_ATOM, s2);
        unsigned nd = rotamer2.get(ROTAMER_ADD_ATOM, s2);

        x[i] = rotamer_angle(atoms[na].X, atoms[nb].X, atoms[nc].X, atoms[nd].X);
      }
      return edge_count;
    }

    template <typename _Atom> unsigned
    write_position_(_I2T<YES_CM_ | YES_UNION_>, real_t *x, const _Atom *atoms,
      const vector_t &angle, const vector_t &ranx) const
    {
      x[0] = angle[0];
      x[1] = angle[1];
      x[2] = angle[2];
      x[3] = ranx[0];
      x[4] = ranx[1];
      x[5] = ranx[2];
      return 6;
    }

  private:

    const _Forcefield *forcefield_; // база данных параметров атомных типов
    const _Residome *residome_;     // база данных топологии
    unsigned freedom_type_;         // предельный уровень свободы

    std::string filename_;          // имя файла загрузки
    std::string molname_;           // имя молекулы (модели)
    unsigned format_;               // формат файла загрузки молекулы
    bool is_solution_;              // признак раствора
    bool exclude_cm_;               // исключать центр инерции при расчетах

    std::vector<_Atomdata>  atomdata_;  // первоначально загруженные атомы
    std::vector<_Bond>      bonds_;     // все связи молекулы
    std::vector<_Angle>     angles_;    // все валентные углы молекулы
    std::vector<_Torsion>   torsions_;  // все дигедралы молекулы
    std::vector<_Pair14>    pair14s_;   // все 1-4 пары молекулы
    std::vector<_Chain >    chains_;    // идентификаторы связных цепей

    std::vector<_Rotamer>   rotamers_;  // идентификаторы атомов, входящих в ротамеры
    std::vector<unsigned>   roots_;     // стартовые ротамеры (их число равно числу несвязанных цепей)
    std::vector<_Edge>      edges_;     // упорядоченные ребра графа ротамеров

    //-----------------------------------------------------------------------------------------------
    // Массив нужнst для внутреннего использования. Их форматы неудобен для внешнего использования,
    // но адаптрованы под эффективную обработку при подготовке данных.
    //-----------------------------------------------------------------------------------------------
    std::set<unsigned>  tor12s_; // номера связей в bonds, которые принадлежат центрам торсионов
    std::set<unsigned>  fix12s_; // номера связей в bonds, которые не вращаются
      // например принадлежат циклам (бензолы, пятичленники) или входят в главную пептидную цепь

    std::vector<unsigned> rotamer_vector_; // хранилище номеров атомов для ротамеров
      // 10, 21, 13, 16, 18, ..., 2219 - номера атомов, связанных с ротамерами в сплошном потоке
      // <... rot1 ...> <... rot2 ...> общее количество элементов равно числу атомов
    std::vector<unsigned> chains_vector_; // хранилище номеров ротамеров для цепей (тот же формат)

    matrix<ACCUMULATE_> t12_matrix_;   // матрица 1-2 связности всех атомов
    matrix<ACCUMULATE_> t13_matrix_;   // матрица 1-3 связности всех атомов
    matrix<ACCUMULATE_> t14_matrix_;   // матрица 1-4 связности всех атомов
      // сохраняется только верхний треугольник матрицы
      // ряды матрицы являются set<int>, чтобы избежать дублирования при построении

      // ВНИМАНИЕ! попытка избежать создания матриц tij_matrix и использовать
      // вместо них сами углы, дигедралы приводит к тому, что часть 1-4 взаимодействий
      // дублируется (например, за счет того что 2 дигедрала могут опираться на
      // одинаковые краевые атомы в ароматических кольцах)

    matrix<ACCUMULATE_> tb_matrix_; // матрица, позволяющая найти все связи атома
      // Она хранит для каждого атома номера связей, в которые он попадает
    matrix<ACCUMULATE_> cr_matrix_; // матрица, позволяющая найти все ротамеры для цепи
      // Она хранит для каждой цепи номера ротамеров, которые в нее попадают
  };

  #define TEMPLATE_HEADER
  #define TEMPLATE_ARG     FORCEFIELD_AMBER_, RESIDOME_AMBER_

  TEMPLATE_HEADER
  INLINE int Archetype_<TEMPLATE_ARG>
  ::find(_I2T<EXT_ID_>, int sid, int pos) const
  {
    int sz = (int)atomdata_.size();
    assert(_GT(sz, 0));

    int lpos = pos, rpos = pos + 1;
    while (lpos >= 0 || rpos < sz)
    {
      if (atomdata_[lpos].sid == sid) return lpos;
      else if (lpos >=0) lpos--;
      if (atomdata_[rpos].sid == sid) return rpos;
      else if (rpos < sz) rpos++;
    }
    return nill;
  }

  TEMPLATE_HEADER
  inline bool Archetype_<TEMPLATE_ARG>
  ::load(const std::string &filename, char altpos)
  {
    filename_ = filename;
    std::string ext = extension(filename);

    std::ifstream file;
    if (ext != _S(".bmm")) file.open(filename.c_str());
    else file.open(filename.c_str(), std::ios_base::binary);
      // bmm файлы, в отличие от других типов файлов, имеют бинарный формат

    if (!file)
    {
      std::string msg = _S("[ERROR] can't open file ") + filename;
      PRINT_BREAK(msg);
    }
    atomdata_.clear();

  #define BUILD_AND_RETURN_TRUE(format) \
    { \
      LOADED_OK_MESSAGE(filename); \
      return true; \
    }

    if ((ext == _S(".pdb") || ext == _S(".ent"))
      && load(_I2T<FORMAT_PDB_>(), file, altpos))
      BUILD_AND_RETURN_TRUE(FORMAT_PDB_)

    if (ext == _S(".hin") && load(_I2T<FORMAT_HIN_>(), file))
      BUILD_AND_RETURN_TRUE(FORMAT_HIN_);

    if (ext == _S(".mol2") && load(_I2T<FORMAT_MOL2_>(), file))
      BUILD_AND_RETURN_TRUE(FORMAT_MOL2_);

    if (ext == _S(".bmm") && load(_I2T<FORMAT_BMM_> (), file))
      BUILD_AND_RETURN_TRUE(FORMAT_BMM_);

  #undef BUILD_AND_RETURN_TRUE

    LOADED_ERR_MESSAGE(filename);
    return false;
  }

  TEMPLATE_HEADER
  INLINE bool Archetype_<TEMPLATE_ARG>
  ::load(_I2T<FORMAT_PDB_>, std::ifstream &file, char altpos)
  {
    molname_ = _S(""); // не знаем имени модели
    format_ = FORMAT_PDB_;
    atomdata_.reserve(DEFAULT_ATOM_COUNT);

    _Atomdata atomdata;
    while (extract_object(_I2T<FORMAT_PDB_>(), atomdata, &file) == CODE_SUCCESS)
    {
      if (atomdata.altloc == ' ' || atomdata.altloc == altpos)
        atomdata_.push_back(atomdata);
    }
    return true;
  }

  TEMPLATE_HEADER
  INLINE bool Archetype_<TEMPLATE_ARG>
  ::load(_I2T<FORMAT_HIN_>, std::ifstream &file, char)
  {
    format_ = FORMAT_HIN_;

    const int DEFAULT_HIN_ATOM_COUNT = 200;
    char skip[256]; skip[0] = 0;
    while (::strncmp(skip, "mol", 3) != 0)
    {
      file.getline(skip, 256);
      if (file.eof()) return false;
    }
    if (skip[3] != '\0') molname_ = _S(&skip[4]);

    atomdata_.reserve(DEFAULT_HIN_ATOM_COUNT);
    _Atomdata atomdata;
    while (extract_object(_I2T<FORMAT_HIN_>(), atomdata, &file) == CODE_SUCCESS)
      atomdata_.push_back(atomdata);

    return true;
  }

  TEMPLATE_HEADER
  INLINE bool Archetype_<TEMPLATE_ARG>
  ::load(_I2T<FORMAT_MOL2_>, std::ifstream &file, char)
  {
    format_ = FORMAT_MOL2_;

    const char *s_molecule = "@<TRIPOS>MOLECULE";
    const char *s_atom = "@<TRIPOS>ATOM";
    const char *s_bond = "@<TRIPOS>BOND";
    char skip[256]; skip[0] = 0;
    char line[256]; line[0] = 0;
    int atom_count = 0, bond_count = 0;
    while (::strncmp(skip, s_molecule, ::strlen(s_molecule)) != 0)
    {
      file.getline(skip, 256);
      if (file.eof()) return false;
    }
    file.getline(skip, 256);
    molname_ = _S(&skip[0]);

    file >> atom_count >> bond_count;
    while (::strncmp(skip, s_atom, ::strlen(s_atom)) != 0)
      file.getline(skip, 256);

    atomdata_.resize(atom_count);
    for (int i=0; i<atom_count; i++)
      extract_object(_I2T<FORMAT_MOL2_>(), atomdata_[i], &file);

    while (::strncmp(skip, s_bond, ::strlen(s_bond)) != 0)
      file.getline(skip, 256);
    int count, atom, atom__; char bondtype__;
    for (int i=0; i<bond_count; i++)
    {
      // sample of format MOL2 line (id, atom, atom, bondtype, options)
      //     1    1    2 2
      file.getline(line, 256);
      ::sscanf(line, "%d %d %d %s", &count, &atom, &atom__, skip);
        // some files after antechamber have rubbish in end of line
        // that produces problems
        // so we treat last field by another way
      skip[2] = '\0';
      if (skip[0] == '1') bondtype__ = 's'; // single
      else if (skip[0] == '2') bondtype__ = 'd'; // double
      else if (skip[0] == '3') bondtype__ = 't'; // triple
      else if (::strncmp(skip, "ar", 2) == 0) bondtype__ = 'a'; // aromatic
      else if (::strncmp(skip, "am", 2) == 0) bondtype__ = 'n'; // amide
      else if (::strncmp(skip, "du", 2) == 0) bondtype__ = 'u'; // dummy
      else if (::strncmp(skip, "nc", 2) == 0) bondtype__ = 'z'; // not connected
      else if (::strncmp(skip, "un", 2) == 0) bondtype__ = '?'; // unknown
      else bondtype__ = '?';

      int ndx = find(_I2T<EXT_ID_>(), atom, atom-1);
      if (ndx == nill)
      {
        _S msg = _S("The file %s has some errors in format", filename_.c_str());
        PRINT_ERR(msg);
      }

      int ndx__ = find(_I2T<EXT_ID_>(), atom__, atom__-1);
      if (ndx__ == nill)
      {
        _S msg = _S("The file %s has some errors in format", filename_.c_str());
        PRINT_ERR(msg);
      }

      _Atomdata &atomdata = atomdata_[ndx];
      _Atomdata &atomdata__ = atomdata_[ndx__];

      atomdata.nid[atomdata.nbond] = atom__;
      atomdata__.nid[atomdata__.nbond] = atom;
      atomdata.nvalency[atomdata.nbond] = bondtype__;
      atomdata__.nvalency[atomdata__.nbond] = bondtype__;
      atomdata.nbond++;
      atomdata__.nbond++;
    }

    return true;
  }

  TEMPLATE_HEADER
  INLINE bool Archetype_<TEMPLATE_ARG>
  ::load(_I2T<FORMAT_BMM_>, std::ifstream &file, char)
  {
    format_ = FORMAT_BMM_;

    file.seekg(sizeof(float) + 3 * sizeof(unsigned));
      // пропустим радиус взаимодействия и размеры ящика

    unsigned atom_count = 0;
    file.read((char*)&atom_count, sizeof(unsigned));
      // прочли число атомов в файле

    atomdata_.reserve(atom_count);
    {
      _Atomdata atomdata;
      while (extract_object(_I2T<FORMAT_BMM_>(), atomdata, &file) == CODE_SUCCESS)
        atomdata_.push_back(atomdata);
    }

    // восстановим обратные ссылки, так как они не храняться в BMM файле
    for (unsigned i=0; i<atom_count; i++)
    {
      _Atomdata &atomdata = atomdata_[i];
      for (unsigned i__=0; i__<atomdata.nbond; i__++)
      {
        if (atomdata.rnid[i__] > 0) // если ссылка вперед
        {
          _Atomdata &atomdata__ = atomdata_[i + atomdata.rnid[i__]];
          atomdata__.rnid[atomdata__.nbond] = -atomdata.rnid[i__];
          atomdata__.nbond++;
        }
      }
    }

    return true;
  }

  TEMPLATE_HEADER
  inline bool Archetype_<TEMPLATE_ARG>
  ::load(_I2T<WATER_>, const std::string &solution)
  {
    filename_ = _S("");
    molname_ = solution; // объявим имя раствора как имя модели
    is_solution_ = true;

    const _Residue *residue = residome_->get_data(solution);
    if (residue)
    {
      typedef _Residue::atom_type _ResidueAtom;
      std::string residue_name = residome_->get_short_name(_I2T<RESIDUE_>(), solution);
      unsigned atom_count = residue->count(_I2T<ATOM_>());

      atomdata_.resize(atom_count);
      for (unsigned i=0; i<atom_count; i++)
      {
        atomdata_[i].sid = i;
        atomdata_[i].residue = residue_name;
      }

      int prev_seq = -1, cur_seq = -1;
        // выделяем cur_seq отдельно, так как в файле могут быть пропуски нумерации
      for (unsigned i=0; i<atom_count; i++)
      {
        _Atomdata &atomdata = atomdata_[i];
        const _ResidueAtom *residue_atom = residue->get_data(i);

        atomdata.pdb_name = residue_atom->name;
        atomdata.fftype = residue_atom->fftype;

        atomdata.charge = residue_atom->charge;
        atomdata.X = residue_atom->X0;

        // вставим связи атомов друг с другом
        atomdata.nbond = residue_atom->na + residue_atom->nh;
        for (unsigned n=0; n<residue_atom->na; n++)
        {
          atomdata.nid[n] = residue_atom->nid[n];
          atomdata.nvalency[n] = 's';
          atomdata.rnid[n] = atomdata.nid[n] - i;
        }
        for (unsigned n=0; n<residue_atom->nh; n++)
        {
          atomdata.nid[n + residue_atom->na] = residue_atom->nhid[n];
          atomdata.nvalency[n + residue_atom->na] = 's';
          atomdata.rnid[n + residue_atom->na] = atomdata.nid[n + residue_atom->na] - i;
        }

        {
          typedef Params<NUCLEAR_, FORCEFIELD_AMBER_> _Param;
          _Param param__;
          _Param::index_type index__(atomdata.fftype);
          forcefield_->get_data(&param__, index__);
          atomdata.mass = param__.mass;

          atomdata.nuclear = find_nuclear(_I2T<MASS_>(), atomdata.mass);
          atomdata.name = nuclears[atomdata.nuclear].name;

        #ifdef USE_HEAVY_HYDROGENS
          if (atomdata.nuclear == 1) atomdata.mass *= HEAVY_MASS_FACTOR;
        #endif
        }

        // построим параметры sigma, так как по ним определяется наличие клешей
        // в процедурах заполнения ящиков водой
        {
          typedef Params<ATOM_, FORCEFIELD_AMBER_> _Param;
          typedef _Param::index_type _Index;
          int nff = forcefield_->get_data_ndx(_I2T<ATOM_>(), _Index(atomdata.fftype));
          if (nff == nill)
          {
            fstring eqv = forcefield_->get(_I2T<EQUI_>(), _Index(atomdata.fftype));
              // get vdw equivalent type
            nff = forcefield_->get_data_ndx(_I2T<ATOM_>(), _Index(eqv));
            if (nff == nill) { NO_EQUIVALENTS(_Index(eqv)); }
          }
          const _Param *param__ = forcefield_->get(_I2T<ATOM_>(), nff);
          atomdata.sigma = param__->sigma;
          atomdata.eps = param__->eps;
        }

        if (residue_atom->resx != prev_seq)
        {
          // при сравнении полагается на то, что в файле _solvent.lib
          // атомы разных цепей не перемещаны друг с другом
          cur_seq++;
          prev_seq = residue_atom->resx;
        }

        atomdata.res_seq = cur_seq; // меняем идентификатор (после cur_seq++)
          // для обеспечения гарантированного упорядочения всех цепей
      }

      LOADED_OK_MESSAGE(solution);
      return true;
    }

    LOADED_ERR_MESSAGE(solution);
    return false;
  }

  TEMPLATE_HEADER
  template <typename _Atom> INLINE unsigned Archetype_<TEMPLATE_ARG>
  ::save(const std::string &filename, const _Atom *atoms, bool prn_hydrogens) const
  {
    std::ofstream file(filename.c_str());
    if (!file)
    {
      std::string msg = _S("[ERROR] can't open file ") + filename;
      PRINT_BREAK(msg);
    }

    unsigned count = 0;
    std::string ext = extension(filename);

    if (ext == _S(".pdb") || ext == _S(".ent"))
      count = save(_I2T<FORMAT_PDB_> (), file, atoms, prn_hydrogens);
    else if (ext == _S(".hin") ) count = save(_I2T<FORMAT_HIN_> (), file, atoms, prn_hydrogens);
    else if (ext == _S(".mol2")) count = save(_I2T<FORMAT_MOL2_>(), file, atoms, prn_hydrogens);
    else if (ext == _S(".bmm") ) count = save(_I2T<FORMAT_BMM_> (), file, atoms, prn_hydrogens);

    if (count) SAVED_OK_MESSAGE(filename);
    return count;
  }

  TEMPLATE_HEADER
  template <typename _Atom>  INLINE unsigned Archetype_<TEMPLATE_ARG>
  ::save(_I2T<FORMAT_PDB_>, std::ofstream &file, const _Atom *atoms,
    bool prn_hydrogens) const
  {
    char line[256];
    file << "MODEL " << molname_ << "\n";

    unsigned atom_count = atomdata_.size();
    for (unsigned i=0; i<atom_count; i++)
    {
      const _Atomdata &atomdata = atomdata_[i];
      if (!prn_hydrogens && atomdata.is_hydrogen()) continue;

      make_string(_I2T<FORMAT_PDB_>(), line, atoms[i], atomdata);
      file << line << std::endl;
    }
    file << "ENDMDL\n";
    return 1;
  }

  TEMPLATE_HEADER
  template <typename _Atom>
  INLINE unsigned Archetype_<TEMPLATE_ARG>
  ::save(_I2T<FORMAT_HIN_>, std::ofstream &file, const _Atom *atoms, bool) const
  {
    char line[256];
    file << "mol " << molname_ << "\n";

    unsigned atom_count = atomdata_.size();
    for (unsigned i=0; i<atom_count; i++)
    {
      make_string(_I2T<FORMAT_HIN_>(), line, atoms[i], atomdata_[i]);
      file << line << std::endl;
    }
    file << "endmol " << molname_ << "\n";
    return 1;
  }

  TEMPLATE_HEADER
  template <typename _Atom>  INLINE unsigned Archetype_<TEMPLATE_ARG>
  ::save(_I2T<FORMAT_MOL2_>, std::ofstream &file, const _Atom *atoms, bool) const
  {
    char line[256];

    unsigned atom_count = atomdata_.size();

    unsigned bond_count = bonds_.size();
    // подсчитаем связи прямо, если код не создавал массива bonds[]
    if (bond_count == 0)
    {
      for (unsigned i=0; i<atom_count; i++) bond_count += atomdata_[i].nbond;
      bond_count >>= 1; // число меньше вдвое, так как сосчитали связи вперед и назад
    }

    file << "@<TRIPOS>MOLECULE\n";
    file << molname_ << '\n';
    ::sprintf(line, "%5d %5d     1     0     0\n", (int)atom_count, (int)bond_count);
    file << line;
    file << "SMALL\n";
    file << "bcc\n\n\n";
    file << "@<TRIPOS>ATOM\n";

    for (unsigned i=0; i<atom_count; i++)
    {
      make_string(_I2T<FORMAT_MOL2_>(), line, atoms[i], atomdata_[i]);
      file << line << std::endl;
    }
    file << "@<TRIPOS>BOND\n";
    int count = 0, nvalency = 0;
    for (unsigned i=0; i<atom_count; i++)
    {
      const _Atomdata &atomdata = atomdata_[i];
      int sid = atomdata.sid;
      for (unsigned i__=0; i__<atomdata.nbond; i__++)
      {
        if (sid > atomdata.nid[i__]) continue; // forward only
        if (atomdata.nvalency[i__] == 's') nvalency = 1;
        if (atomdata.nvalency[i__] == 'd') nvalency = 2;
        if (atomdata.nvalency[i__] == 't') nvalency = 3;
        ::sprintf(line, "%6d %4d %4d %1d", ++count, sid,
          atomdata.nid[i__], nvalency);
        file << line << std::endl;
      }
    }
    return 1;
  }

  TEMPLATE_HEADER
  template <typename _Atom>  INLINE unsigned Archetype_<TEMPLATE_ARG>
  ::save(_I2T<FORMAT_BMM_>, std::ofstream &file, const _Atom *atoms, bool) const
  {
    char line[BMM_RECORD_LEN]; // место для записи атома
    unsigned atom_count = atomdata_.size();
    for (unsigned i=0; i<atom_count; i++)
    {
      make_string(_I2T<FORMAT_BMM_>(), line, atoms[i], atomdata_[i]);
      file.write(line, BMM_RECORD_LEN);
    }
    return atom_count;
  }

  TEMPLATE_HEADER
  inline void Archetype_<TEMPLATE_ARG>
  ::build(_I2T<MOLECULE_>, int pH)
  {
    switch (format_)
    {
    case FORMAT_PDB_  : build_(FORMAT_PDB,  pH); break;
    case FORMAT_HIN_  : build_(FORMAT_HIN,  pH); break;
    case FORMAT_MOL2_ : build_(FORMAT_MOL2, pH); break;
    case FORMAT_BMM_  : build_(FORMAT_BMM,  pH); break;
    }

    unsigned build_components = BUILD_NOTHING;

    if (freedom_type_ & (YES_ATOM_ | YES_ROTAMER_))
    {
      build_components |= BUILD_BONDS;
      build_components |= BUILD_ANGLES;
      build_components |= BUILD_TORSIONS;
      build_components |= BUILD_PAIRS14;
      build_components |= BUILD_CONNECTS;
      build_components |= BUILD_ROTAMERS;
      build_components |= BUILD_CHAINS;
    }

    if (build_components)
    {
      _S msg = _S("Building of model `") + molname_ + _S("` is started ...");
      PRINT_MESSAGE(msg);

      if (build_components & BUILD_BONDS   ) build_(BOND);
      if (build_components & BUILD_ANGLES  ) build_(ANGLE);
      if (build_components & BUILD_TORSIONS) build_(TORSION);
      if (build_components & BUILD_PAIRS14 ) build_(PAIR14);
      if (build_components & BUILD_CONNECTS) build_(CONNECT);
      if (build_components & BUILD_ROTAMERS) build_(ROTAMER);
      if (build_components & BUILD_CHAINS  ) build_(CHAIN);
    }
  }

  TEMPLATE_HEADER
  INLINE bool Archetype_<TEMPLATE_ARG>
  ::build(_I2T<WATER_>, int)
  {
    unsigned build_components = BUILD_NOTHING;

    if (freedom_type_ & (YES_ATOM_ | YES_ROTAMER_))
    {
      build_components |= BUILD_BONDS;
      build_components |= BUILD_ANGLES;
      build_components |= BUILD_TORSIONS;
      build_components |= BUILD_PAIRS14;
      build_components |= BUILD_CONNECTS;
      build_components |= BUILD_ROTAMERS;
      build_components |= BUILD_CHAINS;
    }
    else if (freedom_type_ & YES_CM_)
      build_components |= BUILD_CHAINS;

    if (build_components)
    {
      _S msg = _S("Building of model `") + molname_ + _S("` is started ...");
      PRINT_MESSAGE(msg);

      if (build_components & BUILD_BONDS   ) build_(BOND);
      if (build_components & BUILD_ANGLES  ) build_(ANGLE);
      if (build_components & BUILD_TORSIONS) build_(TORSION);
      if (build_components & BUILD_PAIRS14 ) build_(PAIR14);
      if (build_components & BUILD_CONNECTS) build_(CONNECT);
      if (build_components & BUILD_ROTAMERS) build_(ROTAMER);
      if (build_components & BUILD_CHAINS  ) build_(CHAIN);
    }

    return true;
  }

  TEMPLATE_HEADER
  inline bool Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<FORMAT_PDB_>, int pH)
  {
    build_(RESIDUE_NAME, pH);
    build_(HYDROGEN);

    for (unsigned i=0,sz=atomdata_.size(); i<sz; i++)
      make_object(_I2T<FORMAT_PDB_>(), atomdata_[i], forcefield_, residome_);
      // построение перенесено за восстановлением имен остатков и за
      // построением водородов, чтобы параметры атомов соответствовали
      // именам остатков при pH и вновь вставленные атомы строились
      // в едином месте

    if (freedom_type_ & (YES_ROTAMER_ | YES_ATOM_))
    {
      build_(RESIDUE_CONTACT);
    }
    return true;
  }

  TEMPLATE_HEADER
  inline bool Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<FORMAT_HIN_>, int)
  {
    for (unsigned i=0,sz=atomdata_.size(); i<sz; i++)
      make_object(FORMAT_HIN, atomdata_[i], forcefield_, residome_);

    // создадим относительные ссылки на atomdata индексы
    for (unsigned i=0,sz=atomdata_.size(); i<sz; i++)
    {
      unsigned nbond = atomdata_[i].nbond;
      for (unsigned j=0; j<nbond; j++)
      {
        unsigned nid = atomdata_[i].nid[j];
        int ndx = find(_I2T<EXT_ID_>(), nid, nid-1);
        if (ndx == nill)
        {
          _S msg = _S("The file %s has some errors in format", filename_.c_str());
          PRINT_ERR(msg);
        }
        atomdata_[i].rnid[j] = ndx - i; // ссылки на относительные atomdata индексы
      }
    }

    return true;
  }

  TEMPLATE_HEADER
  INLINE bool Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<FORMAT_MOL2_>, int)
  {
    return build_(FORMAT_HIN);
    // построение полностью идентично
  }

  TEMPLATE_HEADER
  INLINE bool Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<FORMAT_BMM_>, int)
  {
    // строим отсутствующие sid & nid
    for (unsigned i=0,sz=atomdata_.size(); i<sz; i++)
    {
      _Atomdata &atomdata = atomdata_[i];
      atomdata.sid = i + 1;
      for (unsigned i__=0; i__<atomdata.nbond; i__++)
      {
        atomdata.nid[i__] = i + 1 + atomdata.rnid[i__];
        atomdata.nvalency[i__] = 's';
      }
    }

    return build_(FORMAT_HIN);
      // построение полностью идентично
  }

   TEMPLATE_HEADER
   inline void Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<RESIDUE_NAME_>, int pH)
  {
  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_START("Building of residue names : ", 1);
  #endif
    assert(_NE((void*)residome_, (void*)NULL));

    unsigned atom_count = atomdata_.size();
    //--------------------------------------------------------------------------
    //                 корректируем все имена согласно pH
    //--------------------------------------------------------------------------
    std::string residue = "###"; int prev_seq = -1;
    for (unsigned cur=0; cur<atom_count; cur++)
    {
      _Atomdata &atomdata = atomdata_[cur];
      if (atomdata.res_seq != prev_seq)
      {
        residue = residome_->get_equivalent(_I2T<RESIDUE_>(), atomdata.residue, pH);
        prev_seq = atomdata.res_seq;
      }
      atomdata.residue = residue;
    }

    //--------------------------------------------------------------------------
    //     поиск и переименование остатков, которые могут иметь S-S связи
    //--------------------------------------------------------------------------
    typedef std::pair<unsigned, vector_t> S_atom; // <номер атома, координата S>
    std::vector<S_atom> s_atoms;

    for (unsigned cur=0; cur<atom_count; cur++)
    {
      _Atomdata &atomdata = atomdata_[cur];

      // Полагаем, что все атомы серы, находящиеся на расстоянии связи,
      // должны быть связаны друг с другом независимо от имен остатков.
      // Исключением является MET.
      if ( make_string(atomdata.pdb_name).substr(0,1) == _S("S")
        && atomdata.residue != _S("MET"))
      {
        s_atoms.push_back(S_atom(cur, atomdata.X));
      }
    }

    const real_t SS_BOND_LENGTH_2 = sqr(SS_BOND_LENGTH);
    unsigned ssbond_count = 0; // число сульфидных мостиков
    std::vector<unsigned> s_numbers; // номера атомов серы

    // Используем полный перебор пар, так как S-S пар достаточно мало.
    for (unsigned i=0, sz=s_atoms.size(); i<sz; i++)
    {
      vector_t X = s_atoms[i].second;
      for (unsigned j=i+1; j<sz; j++)
      {
        vector_t X__ = s_atoms[j].second;
        if (distance2(X, X__) < SS_BOND_LENGTH_2)
        {
          s_numbers.push_back(s_atoms[i].first);
          s_numbers.push_back(s_atoms[j].first);
          ssbond_count++;
          break; // не может быть нескольких мостиков с одним атомом
        }
      }
    }

    // Заметим, что назначаемые имена "CYX" далее могут быть переименованы,
    // если они попадают в N и C терминалы.
    for (unsigned i=0, sz=s_numbers.size(); i<sz; i++)
    {
      int resid = atomdata_[s_numbers[i]].res_seq;
      int ndx = s_numbers[i];

      std::string res_name = make_string(atomdata_[ndx].residue);
      do atomdata_[ndx--].residue = _S("CYX");
      while (ndx >= 0 && atomdata_[ndx].res_seq == resid);

      ndx = s_numbers[i];
      do atomdata_[ndx++].residue = _S("CYX");
      while (atomdata_[++ndx].res_seq == resid);

      std::string msg = make_string("[WARNING] The name of resigue (%3d %4s) has been changed -> CYX",
        resid, res_name.c_str());
      PRINT_MESSAGE(msg);
    }

    if (ssbond_count)
    {
      std::string msg = make_string("[DEBUG] The model has %d S-S bridges.", ssbond_count);
      PRINT_MESSAGE(msg);
    }

    //--------------------------------------------------------------------------
    //                     находим N & C терминалы
    //--------------------------------------------------------------------------
    std::set<unsigned> term; // список атомов, которые являются первыми для N-term

    char prev_chain = '#'; prev_seq = -2; // гарантируем разрыв (delta > 1)
    for (unsigned cur=0; cur<atom_count; cur++)
    {
      _Atomdata &atomdata = atomdata_[cur];
      if (!residome_->has_residue(atomdata.residue)) continue;
        // игнорируем дополнительные HETATOMS

      // контроль смены цепей
      if (atomdata.chain != prev_chain)
      {
        term.insert(cur);
        prev_chain = atomdata.chain;
      }

      // контроль разрыва цепей
      if (atomdata.res_seq != prev_seq)
      {
        if (abs(atomdata.res_seq - prev_seq) > 1) term.insert(cur);
        prev_seq = atomdata.res_seq;
      }
    }

    //--------------------------------------------------------------------------
    //                   подменяем имена N & C терминалов
    //--------------------------------------------------------------------------
    std::set<unsigned>::iterator sit = term.begin(), site = term.end();
    for (; sit!=site; ++sit)
    {
      unsigned ndx = *sit;

      std::string term_residue = residome_->name(_I2T<NTERM_>(), atomdata_[ndx].residue);
      int resid = atomdata_[ndx].res_seq;
      do { atomdata_[ndx].residue = term_residue; }
      while (ndx < atom_count && atomdata_[++ndx].res_seq == resid);

      ndx = mod((int)(*sit) - 1, (int)atom_count);
      while (!residome_->has_residue(atomdata_[ndx].residue)) ndx = mod((int)(--ndx), (int)atom_count);
        // циклически пройдем до конца со стороны больших номеров и найдем
        // C-term пусть и с другой стороны молекулы (при ndx == 0)

      term_residue = residome_->name(_I2T<CTERM_>(), atomdata_[ndx].residue);
      resid = atomdata_[ndx].res_seq;
      do { atomdata_[ndx].residue = term_residue; }
      while (ndx > (unsigned)0 && atomdata_[--ndx].res_seq == resid);
    }

    {
      std::string msg = make_string("[INFO] The model has %d pairs of N & C-terms",
        term.size());
      PRINT_MESSAGE(msg);
    }

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_FINISH;
  #endif
  }

   TEMPLATE_HEADER
   inline void Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<HYDROGEN_>)
  {
    if (atomdata_.empty()) return;

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_START("Building of hydrogens : ", 1);
  #endif

    std::vector<_Atomdata> ex_atomdata;
      // создаваемая расширенная база, включающая отсутствующие атомы

    const real_t DEGENERACY = 0.1; // контроль линейной зависимости векторов

    std::list<_Atomdata> molecule_residue;
      // список атомов, входящих в аминокислотный остаток молекулы

    typedef std::list<_Atomdata>::iterator _Iterator;
    typedef _Residome::residue_type _Residue;
    typedef _Residue::atom_type _ResidueAtom;
    typedef vector_t::value_type real_type;

    vector_t cterm(0.,0.,0.); // позиция атома из предыдущего остатка, который нужен для
     // того, чтобы правильно разместить водороды у N атома текущего остатка

    std::string prev_residue = atomdata_[0].residue;
    int prev_seq = atomdata_[0].res_seq;
    char prev_chain = atomdata_[0].chain;

    std::set<unsigned> present_atoms_ndx; // присутствующие атомы
    std::set<unsigned> missing_atoms_ndx; // отсутствующие атомы
    std::set<fstring> extra_atoms;   // атомы, не описанные в AMBER *.lib
    std::vector<vector_t> atom_positions; // реальные позиции атомов остатка
    std::vector<vector_t> residome_atom_positions; // позиции атомов остатка

    for (unsigned cur=0,sz=atomdata_.size(); cur<=sz; cur++)
      // используем cur==sz, чтобы обработать последний остаток
    {
      if (cur == sz || atomdata_[cur].res_seq != prev_seq)
      {
        // полагаем, что к этой точке создан список атомов остатка
        // обрабатываем как предыдущий остаток
        if (!residome_->has_residue(prev_residue)) continue;
          // игнорируем дополнительные HETATOMS

        //======================================================================
        //                 обрабатываем текущий остаток
        //======================================================================
        present_atoms_ndx.clear();
        missing_atoms_ndx.clear();
        extra_atoms.clear();

        const _Residue *residome_residue = residome_->get_data(prev_residue);
        assert(_NE((void*)residome_residue, (void*)0));

        unsigned atoms_count = residome_residue->count(_I2T<ATOM_>());
        atom_positions.resize(atoms_count);

        _Iterator it, itb = molecule_residue.begin(), ite = molecule_residue.end();
        for (it=itb; it!=ite; ++it)
        {
          int ndx = residome_residue->get_data_ndx(it->pdb_name);
          if (!it->is_hydrogen() && ndx != nill)
          {
            atom_positions[ndx] = it->X;
            present_atoms_ndx.insert(ndx);
          }
          else
          {
            // удаляем как и атомы, незарегистрированные для данного pH,
            // так и водороды, поскольку их можно восстановить очень точно,
            // и далее они восстанавливаются (так что нужно избежать дублирования
            // при восстановлении
            extra_atoms.insert(it->pdb_name);
            it->pdb_name = "####";
          }
        }

        if (!extra_atoms.empty())
        {
          std::string msg = make_string("[WARNING] Extra atoms of residue %s (%d) [",
            make_string(prev_residue).c_str(), prev_seq);

          std::set<fstring>::iterator it = extra_atoms.begin(),
            ite = extra_atoms.end();
          for (; it!=ite; ++it)
          {
            msg += _S(" ") + make_string(*it);
          }
          msg += _S(" ] have been deleted.");
          PRINT_MESSAGE(msg);

          molecule_residue.remove_if(equal_<_Atomdata, fstring, &_Atomdata::pdb_name>("####"));
            // реальное удаление (так как это список) всех элементов списка, значение поля pdb_name
            // для элементов которого равно "####"
        }

        //----------------------------------------------------------------------
        //               найдем реально отсутствующие атомы
        //----------------------------------------------------------------------
        std::set<unsigned>::iterator sit = present_atoms_ndx.begin(),
          site = present_atoms_ndx.end();
        for (unsigned i=0; i<atoms_count; i++)
        {
          if (sit == site || i != *sit)
          {
            const _ResidueAtom *atom = residome_residue->get_data(i);
            std::string name = make_string(atom->name);

            if (name[0] != 'H') missing_atoms_ndx.insert(i);
              // вставка всех атомов, кроме водородов, для которых
              // координаты вставки могут быть определены значительно точнее
              // вне метода морфинга
            continue;
          }
          if (sit != site && i == *sit) ++sit;
        }

        _Atomdata atomdata__;
        atomdata__.residue = prev_residue;
        atomdata__.chain = prev_chain;
        atomdata__.res_seq = prev_seq;

        if (!missing_atoms_ndx.empty())
        {
          std::string msg = make_string("[WARNING] The residue %s (%d %c) have next "
          "missing atoms to be restored : ", prev_residue.c_str(), prev_seq, prev_chain);

          std::set<unsigned>::iterator it = missing_atoms_ndx.begin(),
            ite = missing_atoms_ndx.end();
          for (; it!=ite; ++it)
          {
            const _ResidueAtom *atom = residome_residue->get_data(*it);
            msg += make_string(atom->name) + _S(" ");
          }
          PRINT_MESSAGE(msg);

          //--------------------------------------------------------------------
          //                   восстановим пропущенный OXT
          //--------------------------------------------------------------------
            //  Данное восстановление наиболее часто требуется, по этому выделено
            //  как отдельное с большей точностью восстановления.

          if (missing_atoms_ndx.size() == 1
            && make_string(residome_residue->get_data(*missing_atoms_ndx.begin())->name) == _S("OXT")
          )
          {
            vector_t XA[4] = {0., 0., 0., 0.}, XS[4];
            XS[0] = residome_residue->get(_I2T<ATOM_>(), "C"  )->X0;
            XS[1] = residome_residue->get(_I2T<ATOM_>(), "CA" )->X0;
            XS[2] = residome_residue->get(_I2T<ATOM_>(), "O"  )->X0;
            XS[3] = residome_residue->get(_I2T<ATOM_>(), "OXT")->X0;

            _Iterator it = molecule_residue.begin(), ite = molecule_residue.end();
            for (; it!=ite; ++it)
            {
              const _Atomdata &atomdata = *it;
              _S pdb_name = make_string(atomdata.pdb_name);
              if (pdb_name == _S("C" )) XA[0] = atomdata.X;
              if (pdb_name == _S("CA")) XA[1] = atomdata.X;
              if (pdb_name == _S("O" )) XA[2] = atomdata.X;
            }

            build_2A_1A_(&XA[0], &XS[0]);

            // сохранение пропущенных атомов
            //------------------------------
            atomdata__.X = XA[3];
            atomdata__.pdb_name = "OXT";
            molecule_residue.push_back(atomdata__);
          }
          else
          {
            //----------------------------------------------------------------------
            //                восстановим другие пропущенные атомы
            //----------------------------------------------------------------------
            unsigned n = present_atoms_ndx.size();
            if (n < 3)
            {
              std::string msg = make_string("[ERROR] It's impossible to restore atoms. "
                "Residue %s (%d) has a few atoms.\n", prev_residue.c_str(), prev_seq);
              PRINT_BREAK(msg);
            }

            residome_atom_positions.resize(atoms_count);
            for (unsigned i=0; i<atoms_count; i++)
            {
              const _ResidueAtom *atom = residome_residue->get_data(i);
              residome_atom_positions[i] = atom->X0;
            }

            typedef std::pair<unsigned, real_t> _Pair;
            std::vector<_Pair> R2(n);
              // { index of present atom, distance2 between MISSING & PRESENT atoms }

            std::set<unsigned>::iterator mit, pit,
              mitb = missing_atoms_ndx.begin(),
              mite = missing_atoms_ndx.end(),
              pitb = present_atoms_ndx.begin(),
              pite = present_atoms_ndx.end();

            for (mit=mitb; mit!=mite; ++mit)
            {
              vector_t X0 = residome_atom_positions[*mit];
              pit = pitb;
              for (unsigned i=0; pit!=pite; ++i, ++pit)
              {
                vector_t X = residome_atom_positions[*pit];
                R2[i].first  = *pit;
                R2[i].second = distance2(X, X0);
              }

              std::sort(R2.begin(), R2.end(), less_<_Pair, real_t, &_Pair::second>());
                // рассчитываем и упорядочиваем все расстояния между
                // MISSING ATOM & PRESENT ATOMS

            #ifdef FULL_BUILD_DEBUG
              {
                const _ResidueAtom *atom = residome_residue->get_data(*mit);
                std::string msg = _S("[DEBUG] ") + make_string(prev_residue)
                  + _S("(") + itoa(prev_seq) + _S(") has next sort based atoms for missing : ")
                  + make_string(atom->name) + _S("\n");

                for (unsigned i=0, sz=R2.size(); i<sz; i++)
                {
                  unsigned ndx = R2[i].first;
                  const _ResidueAtom *atom__ = residome_residue->get_data(ndx);
                  msg += _S("  ") + make_string(atom__->name)
                    + _S("/") + ftoa(sqrt(R2[i].second)) + _S("/");
                }
                PRINT_MESSAGE(msg);
              }
            #endif

              // выбор базисных векторов морфинга
              //---------------------------------

              unsigned nn = 0; // индексы в сортированном r2 массиве присутствующих атомов
              vector_t A, B, C, // локальная аминокислотная система
                A__, B__, C__; // молекулярная система
              int ndxA, ndxB, ndxC; // индексы атомов в локальной системе
              real_t sp = 0;

              ndxA = R2[nn++].first;
              A = residome_atom_positions[ndxA];

              ndxB = R2[nn++].first;
              B = residome_atom_positions[ndxB];

              do {
                ndxC = R2[nn++].first;
                C = residome_atom_positions[ndxC];
                vector_t AB = A - B;
                vector_t AC = A - C;
                sp = sqr(scalar_product(AB, AC)) - AB.length2() * AC.length2();
                  // контроль, чтобы атомы не оказались на одной прямой
              }
              while (fabs(sp) < DEGENERACY && nn < n);

              if (nn >= n)
              {
                PRINT_ERR("Can't restore missing atoms. Degeneracy of local system.\n");
              }

              A__ = atom_positions[ndxA];
              B__ = atom_positions[ndxB];
              C__ = atom_positions[ndxC];

              // Делаем совмещение аминокислотной системы в молекулярную по трем точкам.
              // В результате образ атома аминокислотной системы (отсутствующий в молекуле)
              // попадет на предполагаемое место, где должен быть атом в молекуле.

              // совмещаем молекулярную и аминокислотную координатные системы по точке A
              vector_t AA = A__ - A; X0 += AA; A = A__; B += AA; C += AA;

              // совмещаем молекулярную и аминокислотную координатные системы по оси AB
              {
                vector_t BA = B - A__, BA__ = B__ - A__;
                vector_t W = vector_product(BA, BA__); // ось вращения
                real_t w = get_angle(BA, BA__); // угод вращения для аминокислотной системы

                Rotator<AXIS_ROTATOR_, real_t> rotator(W, w, A__);
                X0 = rotator(X0); A = rotator(A); B = rotator(B); C = rotator(C);
              }

              // совмещаем молекулярную и аминокислотную координатные системы по плоскости ABC
              {
                vector_t W = B - A__; W.normalize(); // ось вращения
                vector_t CA = C - A__, CA__ = C__ - A__;
                vector_t CP = CA  - W * scalar_product(W, CA),
                  CP__ = CA__ - W * scalar_product(W, CA__);
                  // перпендикуляры к оси вращения W

                real_t w = get_angle(CP, CP__); // угод вращения для аминокислотной системы

                Rotator<AXIS_ROTATOR_, real_t> rotator(W, w, A__);
                X0 = rotator(X0); A = rotator(A); B = rotator(B); C = rotator(C);
              }

              atom_positions[*mit] = X0;

            #ifdef FULL_BUILD_DEBUG
              {
                const _ResidueAtom *atom__ = residome_residue->get_data(*mit);
                std::string msg = _S("  position of missing atom ")
                  + make_string(atom__->name) + _S(" ") + make_string(atom_positions[*mit]);
                PRINT_MESSAGE(msg);
              }
            #endif
            }

            // сохранение пропущенных атомов
            //------------------------------
            for (mit=mitb; mit!=mite; ++mit)
            {
              const _ResidueAtom *atom__ = residome_residue->get_data(*mit);

              atomdata__.X = atom_positions[*mit];
              atomdata__.pdb_name = atom__->name;

              molecule_residue.push_back(atomdata__);
            }
          }
        }
        //----------------------------------------------------------------------
        //                  конец восстановления атомов
        //----------------------------------------------------------------------


        //----------------------------------------------------------------------
        //                    восстановление водородов
        //----------------------------------------------------------------------
        // предполагаем, что атомов водорода нет в исходном PDB файле
        vector_t X0; vector_t XH[3], XA[3];

        for (it=itb; it!=ite; ++it)
        {
          int ndx0 = residome_residue->get_data_ndx((*it).pdb_name);
          const _ResidueAtom *atom = residome_residue->get_data(ndx0);

          unsigned sza = atom->na, szh = atom->nh;
          if (szh == 0) continue; // пропустим атом без водородов

          X0 = (*it).X;

          // рассмотрим соседей атома, к кому присоединяются водороды
          for (unsigned i=0; i<sza; i++)
          {
            const _ResidueAtom *atom__ = residome_residue->get_data(atom->nid[i]);
            for (_Iterator it__=itb; it__!=ite; ++it__)
            {
              if ((*it__).pdb_name == atom__->name)
              { XA[i] = (*it__).X; break; }
            }
          }
          if (residome_residue->get(_I2T<EXT_CONTACT_PREV_>()) == ndx0)
            XA[sza++] = cterm; // требуется внешний атом

          if (sza == 1)
          {
            if (szh == 1) build_1A_1H_(&XH[0], &XA[0], X0, (real_type)HYDROGEN_ATOM_DISTANCE);
            if (szh == 2) build_1A_2H_(&XH[0], &XA[0], X0, (real_type)HYDROGEN_ATOM_DISTANCE);
            if (szh == 3) build_1A_3H_(&XH[0], &XA[0], X0, (real_type)HYDROGEN_ATOM_DISTANCE);
          }
          if (sza == 2)
          {
            if (szh == 1) build_2A_1H_(&XH[0], &XA[0], X0, (real_type)HYDROGEN_ATOM_DISTANCE);
            if (szh == 2) build_2A_2H_(&XH[0], &XA[0], X0, (real_type)HYDROGEN_ATOM_DISTANCE);
          }
          if (sza == 3 && szh == 1)
            build_3A_1H_(&XH[0], &XA[0], X0, (real_type)HYDROGEN_ATOM_DISTANCE);

          for (unsigned i=0; i<szh; i++)
          {
            const _ResidueAtom *atom__ = residome_residue->get_data(atom->nhid[i]);
            atomdata__.X = XH[i];
            atomdata__.pdb_name = atom__->name;

            molecule_residue.push_back(atomdata__);
          }
        } // конец восстановления водородов

        //======================================================================
        //             переходим к обработке следующего остатка
        //======================================================================
        if (cur < sz)
        {
          prev_residue = atomdata_[cur].residue;
          prev_seq = atomdata_[cur].res_seq;
          prev_chain = atomdata_[cur].chain;

          // восстановим позицию атома, необходимого для связи со следующим остатком
          int contact = residome_residue->get(_I2T<EXT_CONTACT_NEXT_>());
          if (contact != nill)
          {
            const _ResidueAtom *atom__ = residome_residue->get_data(contact);
            _Iterator it = itb;
            for (; it!=ite; ++it)
              if ((*it).pdb_name == atom__->name) { cterm = (*it).X; break; }
            if (it == ite)
            {
              _S msg = _S("Can't find contact from previous residue.");
              PRINT_ERR(msg);
            }
          }
        }

        // переупорядочим атомы остатка согласно порядку в базе данных residome,
        // что обеспечивает нам минимальность смещений номеров для связанных
        // друг с другом атомов, сохраним остаток во временной базе и очистим его

        _Iterator pos = molecule_residue.begin(); // итератор перед которым вставляю
          // элементы после их нахождения, формируя таким образом сортированный список
        for (unsigned i=0; i<atoms_count; i++)
        {
          const _ResidueAtom *atom__ = residome_residue->get_data(i);
          if ((*pos).pdb_name == atom__->name) { ++pos; continue; }
            // итератор перед которым вставляем и нужно вставлять,
            // просто сдвинемся к следующему без вставки

          for (it=pos; it!=molecule_residue.end(); ++it)
            // конец после вычленения элемента может измениться, потому
            // все время перевычисляется
          {
            if ((*it).pdb_name == atom__->name)
              molecule_residue.splice(pos, molecule_residue, it);
          }
        }

        it = molecule_residue.begin(); ite = molecule_residue.end();
        for (; it!=ite; ++it) ex_atomdata.push_back(*it);
        molecule_residue.clear();
      }

      if (cur < sz) molecule_residue.push_back(atomdata_[cur]);
    }
    atomdata_.swap(ex_atomdata); // сохранение новых данных в базе путем обмена

    {
      std::string msg = make_string("[WARNING] The protein has %d atoms.", atomdata_.size());
      PRINT_MESSAGE(msg);
    }
  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_FINISH;
  #endif
  }

  TEMPLATE_HEADER
  inline void Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<RESIDUE_CONTACT_>)
  {
    assert(_NE((void*)residome_, (void*)NULL));
    if (atomdata_.size() == 0) return;

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_START("Building of residue contacts : ", 1);
  #endif
    std::string residue;
    int prev_seq = nill;

    for (int i=0,sz=(int)atomdata_.size(); i<sz; i++)
      atomdata_[i].nbond = 0; // очистка старых контактов

    //--------------------------------------------------------------------------
    //                        build residue contacts
    //--------------------------------------------------------------------------
    typedef _Residue::atom_type _ResidueAtom;

    const _Residue *cur_residue = NULL;
    const _ResidueAtom *cur_residue_atom = NULL;
    std::string prev_residue = _S("####");
    fstring contact_atom = "####";
    int contact_id = 0, na = 0;
    int contacts[MAX_ATOM_BOND]; // all contacts of atom (atoms + hydrogens)
    prev_seq = nill;

    for (int i=0,sz=(int)atomdata_.size(); i<sz; i++)
    {
      _Atomdata &atomdata = atomdata_[i];
      if (atomdata_[i].res_seq != prev_seq)
      {
        cur_residue = residome_->get_data(atomdata.residue);
        assert(_NE((void*)cur_residue, (void*)0));

        contact_id = cur_residue->get(_I2T<EXT_CONTACT_NEXT_>());
        if (contact_id >= 0)
        {
          cur_residue_atom = cur_residue->get_data(contact_id);
          contact_atom = cur_residue_atom->name;
        }
        prev_seq = atomdata_[i].res_seq;
      }

      na = cur_residue->get_data_ndx(atomdata.pdb_name);
      assert(_NE(na, nill));

      cur_residue_atom = cur_residue->get_data(na);

      //------------------------------------------------------------------------
      //              делаем локальные контакты текущего атома
      //------------------------------------------------------------------------
      int contact_count = cur_residue_atom->get_contacts((int*)&contacts);
      for (int i__=0; i__<contact_count; i__++)
      {
        const _ResidueAtom *atom__ = cur_residue->get_data(contacts[i__]);

        char valency; int c1 = cur_residue->get(_I2T<EXT_CONTACT_PREV_>()), c2 = contact_id;
        // установим валентность в зависимости от того, находятся ли атомы в пептидной цепи
        // установка высокой жесткости позволяет далее не рассматривать данные связи как ротамеры
        // тем самым выделить отдельно пептидную цепь (концы NH3, COO- являются ротамерами)
        if ( (atomdata.pdb_name == "N"  && c1 >= 0 && atom__->name == "CA")
          || (atomdata.pdb_name == "C"  && c2 >= 0 && atom__->name == "CA")
          || (atomdata.pdb_name == "CA" && atom__->name == "N" && c1 >= 0 )
          || (atomdata.pdb_name == "CA" && atom__->name == "C" && c2 >= 0 )
        ) valency = 'd';
        else valency = 's';

        int k = i + 1; // старт поиска вперед
        while (k < sz && atomdata_[k].pdb_name != atom__->name
          && atomdata_[k].res_seq == atomdata.res_seq) k++;
        // если найден атом впереди, записываем данные
        if (k != sz && atomdata_[k].res_seq == atomdata.res_seq)
        {
          atomdata.nid[atomdata.nbond] = atomdata_[k].sid;
          atomdata.rnid[atomdata.nbond] = k - i;
          atomdata.nvalency[atomdata.nbond] = valency;
          atomdata.nbond++;
          atomdata_[k].nid[atomdata_[k].nbond] = atomdata.sid;
          atomdata_[k].rnid[atomdata_[k].nbond] = i - k;
          atomdata_[k].nvalency[atomdata_[k].nbond] = valency;
          atomdata_[k].nbond++;
          continue; // возможно несколько контактов
        }

        k = i - 1; // старт поиска назад для контроля
          // если находим, значит все нормально и ничего не делаем, так как
          // все сделано при поиске вперед, если нет, печатаем сообщение
        while (k >= 0 && atomdata_[k].pdb_name != atom__->name
          && atomdata_[k].res_seq == atomdata.res_seq) k--;
        if (k >= 0 && atomdata_[k].res_seq == atomdata.res_seq)
          continue; // найден

        std::string msg = _S("can't find the internal residue contact\n")
          + _S("  ") + make_string(atomdata) + _S("  <--->  ")
          + make_string(atom__->name) + _S("\n");
        PRINT_ERR(msg);
      }

      //------------------------------------------------------------------------
      //            делаем межаминокислотный контакт текущего атома
      //------------------------------------------------------------------------
      if (atomdata.pdb_name == contact_atom && contact_id >= 0)
      {
        fstring contact_atom__ = "####";
        int k = i + 1; // старт поиска вперед
        while (k < sz && atomdata_[k].res_seq == atomdata.res_seq) k++;
        if (k == sz)
        {
          std::string msg = _S("\ncan't find the next residue\n")
            + _S("  ") + make_string(atomdata);
          PRINT_ERR(msg);
        }

        const _Residue *cur_residue__ = residome_->get_data(atomdata_[k].residue);
        assert(_NE((void*)cur_residue__, (void*)0));

        int contact_id__ = cur_residue__->get(_I2T<EXT_CONTACT_PREV_>());
        assert(_GE(contact_id__, 0));

        const _ResidueAtom *atom__ = cur_residue__->get_data(contact_id__);
        while (k < sz && atomdata_[k].pdb_name != atom__->name) k++;
          // пропускаем неконтактные атомы
        if (k < sz && atomdata_[k].pdb_name == atom__->name)
        {
          atomdata.nid[atomdata.nbond] = atomdata_[k].sid;
          atomdata.rnid[atomdata.nbond] = k - i;
          atomdata.nvalency[atomdata.nbond] = 'd'; // жесткие межаминокислотные связи
          atomdata.nbond++;
          atomdata_[k].nid[atomdata_[k].nbond] = atomdata.sid;
          atomdata_[k].rnid[atomdata_[k].nbond] = i - k;
          atomdata_[k].nvalency[atomdata_[k].nbond] = 'd';
          atomdata_[k].nbond++;
        }
        else
        {
          std::string msg = _S("can't find the external residue contact\n")
            + make_string(atomdata) + _S("\n")
            + make_string(atomdata_[k]);
          PRINT_ERR(msg);
        }
      }
    #ifdef FULL_TOPOLOGY_DEBUG
      {
        std::string msg = make_string(atomdata.name) + _S(" ")
          + make_string(atomdata);
        PRINT_MESSAGE(msg);
      }
    #endif
    }

    //--------------------------------------------------------------------------
    //                    восстановление S-S связей
    //--------------------------------------------------------------------------
    std::vector<unsigned> s_atoms; // номер атома серы

    // найдем все остатки с атомами серы
    for (unsigned i=0,sz=atomdata_.size(); i<sz; i++)
    {
      std::string residue_name = make_string(atomdata_[i].residue);
      if ( make_string(atomdata_[i].pdb_name).substr(0, 1) == _S("S")
        && (residue_name == _S("CYX") || residue_name == _S("NCYX") || residue_name == _S("CCYX"))
      )
      s_atoms.push_back(i);
    }

    // свяжем все остатки с атомами серы (при малых расстояниях S-S)
    const real_t SS_BOND_LENGTH_2 = sqr(SS_BOND_LENGTH);
    for (unsigned i=0, sz=s_atoms.size(); i<sz; i++)
    {
      _Atomdata &atomdata = atomdata_[s_atoms[i]];
      vector_t &X = atomdata.X;
      for (unsigned j=i+1; j<sz; j++)
      {
        _Atomdata &atomdata__ = atomdata_[s_atoms[j]];
        vector_t &X__ = atomdata__.X;
        if (distance2(X, X__) < SS_BOND_LENGTH_2)
        {
          atomdata.nid[atomdata.nbond] = atomdata_[s_atoms[j]].sid;
          atomdata.rnid[atomdata.nbond] = s_atoms[j] - s_atoms[i];
          atomdata.nvalency[atomdata.nbond] = 's';
          atomdata.nbond++;
          atomdata__.nid[atomdata__.nbond] = atomdata_[s_atoms[i]].sid;
          atomdata__.rnid[atomdata__.nbond] = s_atoms[i] - s_atoms[j];
          atomdata__.nvalency[atomdata__.nbond] = 's';
          atomdata__.nbond++;
        }
      }
    }

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_FINISH;
  #endif
  }

  TEMPLATE_HEADER
  inline void Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<BOND_>)
  {
    assert(_NE((void*)forcefield_, (void*)0));

    int sz = (int)atomdata_.size();
    if (sz < 2) return;
      // невозможно построить углы для молекулы из двух атомов

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_START("Building of bonds : ", 1);
  #endif

    t12_matrix_.resize(sz); // build t12_matrix also
    tb_matrix_.resize(sz); // build tb_matrix also

    // preliminary calculation of bonds count (to save time)
    int bond_count = 0;
    {
      for (int i=0; i<sz; ++i)
      {
        _Atomdata &atomdata = atomdata_[i];
        bond_count += atomdata.nbond;
      }
    }
    bond_count >>= 1; // number of bonds

    bonds_.resize(bond_count);
    int ibond = 0;
    for (int i=0; i<sz; ++i)
    {
      const _Atomdata &atomdata = atomdata_[i];
      for (unsigned n=0; n<atomdata.nbond; n++)
      {
        int i__ = i + atomdata.rnid[n];
        if (i__ < i) continue; // игнорируем обратные ссылки

        if (valency_order(atomdata.nvalency[n]) > 1.f)
        {
          // фиксируем связь как жесткую, если ее порядок высокий
          fix12s_.insert(ibond);
        }

        _Bond &bond = bonds_[ibond];
        bond.ndx = index_<2>(i, i__);
        ffindex2 index(atomdata.fftype, atomdata_[i__].fftype);

        if (!forcefield_->get_data(&bond.param, index))
        {
        #ifdef FULL_TOPOLOGY_DEBUG
          int extid   = get(_I2T<EXT_ID_>(), i  );
          int extid__ = get(_I2T<EXT_ID_>(), i__);
          std::string msg = _S("[WARNING] Missing bond ")
            + make_string(index[0]) + _S("[") + _S(itoa(extid  )) + _S("] - ")
            + make_string(index[1]) + _S("[") + _S(itoa(extid__)) + _S("]");
          PRINT_MESSAGE(msg);
        #endif
          if (!forcefield_->get_equivalent_data(&bond.param, index))
          {
            NO_EQUIVALENTS(index);
            throw std::exception();
          }
        }

        // build t12_matrix
        {
          index_<2> ndx(i, i__);
          if (ndx[0] > ndx[1]) swap(ndx);
          t12_matrix_[ndx[0]].insert(ndx[1]);
        }

        // build tb_matrix
        tb_matrix_[i  ].insert(ibond);
        tb_matrix_[i__].insert(ibond);
        ibond++;
      }
    }
    PRINT_MESSAGE(_S("Have been built -> bonds (") + itoa(bonds_.size()) + _S(")"));

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_FINISH;
  #endif
  }

  TEMPLATE_HEADER
  inline void Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<ANGLE_>)
  {
    unsigned atom_count = atomdata_.size();
    if (atom_count < 3) return;
      // невозможно построить углы для молекулы из двух атомов

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_START("Building of angles : ", 1);
  #endif

    t13_matrix_.resize(atom_count);

    typedef matrix<ACCUMULATE_>::row_iterator _Iterator;

    unsigned bond_count = bonds_.size();
    unsigned angle_count = 0;
    for (unsigned i=0; i<bond_count; i++)
    {
      _Bond &bond = bonds_[i];
      for (unsigned I=0; I<2; I++)
      {
        unsigned n1 = bond.ndx[I];
        for (_Iterator it=tb_matrix_[n1].begin(),
          ite=tb_matrix_[n1].end(); it!=ite; ++it)
        {
          if ((unsigned)(*it) <= i) continue; // skip the same & backward bonds
          angle_count++;
        }
      }
    }

    angles_.resize(angle_count);
    int iangle = 0; // current (flick between langle & rangle)
    for (unsigned i=0; i<bond_count; i++)
    {
      _Bond &b1 = bonds_[i];

      for (unsigned I=0; I<2; I++)
      {
        unsigned n1 = b1.ndx[    I];
        unsigned n2 = b1.ndx[1 - I];
        fstring a1 = atomdata_[n1].fftype;
        fstring a2 = atomdata_[n2].fftype;

        for (_Iterator it=tb_matrix_[n2].begin(),
          ite=tb_matrix_[n2].end(); it!=ite; ++it)
        {
          if ((unsigned)(*it) <= i) continue; // skip the same & backward bonds

          _Bond &b2 = bonds_[*it];
          unsigned n3 = b2.ndx[0];
          if (n3 == n2) n3 = b2.ndx[1];
          fstring a3 = atomdata_[n3].fftype;

          _Angle &vangle = angles_[iangle];
          vangle.ndx[0] = n1;
          vangle.ndx[1] = n2;
          vangle.ndx[2] = n3;

          ffindex3 index(a1, a2, a3);
          if (!forcefield_->get_data(&vangle.param, index))
          {
          #ifdef FULL_TOPOLOGY_DEBUG
            int extidA = get(_I2T<EXT_ID_>(), n1);
            int extidB = get(_I2T<EXT_ID_>(), n2);
            int extidC = get(_I2T<EXT_ID_>(), n3);
            std::string msg = _S("[WARNING] Missing angle ")
              + make_string(index[0]) + _S("[") + _S(itoa(extidA)) + _S("] - ")
              + make_string(index[1]) + _S("[") + _S(itoa(extidB)) + _S("] - ")
              + make_string(index[2]) + _S("[") + _S(itoa(extidC)) + _S("]");
            PRINT_MESSAGE(msg);
          #endif
            if (!forcefield_->get_equivalent_data(&vangle.param, index))
            {
              NO_EQUIVALENTS(index);
              throw std::exception();
            }
          }

          {
            index_<2> ndx(n1, n3);
            if (ndx[0] > ndx[1]) swap(ndx);
            t13_matrix_[ndx[0]].insert(ndx[1]);
          }
          iangle++;
        }
      }
    }
    PRINT_MESSAGE(_S("Have been built -> angles (") + itoa(angles_.size()) + _S(")"));

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_FINISH;
  #endif
  }

  TEMPLATE_HEADER
  inline void Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<TORSION_>)
  {
    unsigned atom_count = atomdata_.size();
    t14_matrix_.resize(atom_count);

    if (atom_count < 4) return;
      // невозможно построить углы для молекулы из трех атомов

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_START("Building of torsions : ", 1);
  #endif

    fstring a1, a2, a3, a4; // идентификаторы типов атомов
    unsigned t1, t2, t3, t4; // номера атомов

    typedef matrix<ACCUMULATE_>::row_iterator _Iterator;
    typedef std::pair<unsigned, unsigned> _Index;

    std::set<_Index> bad14_matrix__;
      // пары атомов, на которые опирается более чем один торсион

    unsigned bond_count = bonds_.size();
    int tors_count = 0;
    for (unsigned i=0; i<bond_count; i++)
    {
      _Bond &b2 = bonds_[i];
      t2 = b2.ndx[0];
      t3 = b2.ndx[1];

      _Iterator it = tb_matrix_[t2].begin(), ite = tb_matrix_[t2].end();
      for (; it!=ite; ++it)
      {
        if ((unsigned)(*it) == i) continue; // исключим текущую связь

        _Bond &b1 = bonds_[*it];
        t1 = b1.ndx[0];
        if (t1 == t2) t1 = b1.ndx[1]; // выберем отдаленный атом

        _Iterator it__ = tb_matrix_[t3].begin(), ite__ = tb_matrix_[t3].end();
        for (; it__!=ite__; ++it__)
        {
          if ((unsigned)(*it__) == i) continue; // исключим текущую связь

          _Bond &b2__ = bonds_[*it__];
          t4 = b2__.ndx[0];
          if (t4 == t3) t4 = b2__.ndx[1];

          if (t4 == t1) continue; // совпадение краевых атомов торсиона
            // такое может быть, например для SCP модели воды

          // Строим матрицу 1-4 взаимодействий и записываем пары атомов,
          // на которые опираются больше чем один торсион (для их исключения в дальнейшем)
          unsigned t1__ = t1, t4__ = t4;
          if (t1__ > t4__) std::swap(t1__, t4__);

          std::pair<_Iterator, bool> result = t14_matrix_[t1__].insert(t4__);
          if (result.second == true)
          {
            tors_count++;
              // Здесь происходит излишнее суммирование, если на данную пару атомов
              // далее булет опираться не один торсион. Потому в конце цикла излишние суммы
              // удаляются.
          }
          if (result.second == false || t13_matrix_[t1__].count(t4__))
          {
            bad14_matrix__.insert(_Index(t1__, t4__));
            tors_count--;
              // Здесь уменьшаем ранее суммированное число торсионов (tors_count++).
              // Такой фокус возможен, так как число торсионов, опирающихся на одну пару атомов
              // не может быть более двух. Иначе такое уменьшение нужно делать после цикла.
          }
        }
      }
    }

    if (tors_count)
    {
      torsions_.resize(tors_count);

      int idihe = 0;

      for (unsigned i=0; i<bond_count; i++)
      {
        _Bond &b2 = bonds_[i];

        t2 = b2.ndx[0];
        t3 = b2.ndx[1];
        a2 = atomdata_[t2].fftype;
        a3 = atomdata_[t3].fftype;

        _Iterator it = tb_matrix_[t2].begin(), ite = tb_matrix_[t2].end();
        for (; it!=ite; ++it)
        {
          if ((unsigned)(*it) == i) continue;

          _Bond &b1 = bonds_[*it];
          t1 = b1.ndx[0];
          if (t1 == t2) t1 = b1.ndx[1];
            // выберем отдаленный атом

          a1 = atomdata_[t1].fftype;

          _Iterator it__ = tb_matrix_[t3].begin(), ite__ = tb_matrix_[t3].end();
          for (; it__!=ite__; ++it__)
          {
            if ((unsigned)(*it__) == i) continue;

            _Bond &b2__ = bonds_[*it__];
            t4 = b2__.ndx[0];
            if (t4 == t3) t4 = b2__.ndx[1];

            unsigned t1__ = t1, t4__ = t4;
            if (t1__ > t4__) std::swap(t1__, t4__);
            if (bad14_matrix__.count(_Index(t1__, t4__)))
            {
              fix12s_.insert(i); // записываем номера связей, попавших в бензольные кольца
              continue;
            }

            a4 = atomdata_[t4].fftype;

            _Torsion &tors = torsions_[idihe];
            tors.ndx[0] = t1;
            tors.ndx[1] = t2;
            tors.ndx[2] = t3;
            tors.ndx[3] = t4;

            ffindex4 index(a1, a2, a3, a4);
            if (!forcefield_->get_data(&tors.param, index))
            {
            #ifdef FULL_TOPOLOGY_DEBUG
              int extid1 = get(_I2T<EXT_ID_>(), t1);
              int extid2 = get(_I2T<EXT_ID_>(), t2);
              int extid3 = get(_I2T<EXT_ID_>(), t3);
              int extid4 = get(_I2T<EXT_ID_>(), t4);
              std::string msg = _S("[WARNING] Missing torsion ")
                + make_string(index[0]) + _S("[") + _S(itoa(extid1)) + _S("] - ")
                + make_string(index[1]) + _S("[") + _S(itoa(extid2)) + _S("] - ")
                + make_string(index[2]) + _S("[") + _S(itoa(extid3)) + _S("] - ")
                + make_string(index[3]) + _S("[") + _S(itoa(extid4)) + _S("]");
              PRINT_MESSAGE(msg);
            #endif
              if (!forcefield_->get_equivalent_data(&tors.param, index))
              {
                NO_EQUIVALENTS(index);
                throw std::exception();
              }
            }

            tor12s_.insert(i); // записываем номера связей, попавших в центр торсиона
            idihe++;
          }
        }
      }
    }
    PRINT_MESSAGE(_S("Have been built -> torsions (") + itoa(torsions_.size()) + _S(")"));

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_FINISH;
  #endif
  }

  TEMPLATE_HEADER
  inline void Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<PAIR14_>)
  {
    unsigned atom_count = atomdata_.size();

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_START("Building of pair14s : ", 1);
  #endif

    // сосчитаем число элементов, чтобы знать количество требуемой памяти
    unsigned count = 0;
    for (unsigned i=0; i<atom_count; i++) count += t14_matrix_[i].size();
    pair14s_.resize(count);

    // подготовим массив 1-4 пар
    _Pair14 pair__;
    for (unsigned i=0,k=0; i<atom_count; i++)
    {
      pair__[0] = i;
      std::set<int>::iterator it = t14_matrix_[i].begin(), ite = t14_matrix_[i].end();
      for (; it!=ite; ++it)
      {
        pair__[1] = *it;
        pair14s_[k++] = pair__;
      }
    }

    PRINT_MESSAGE(_S("Have been built -> 1-4 pairs (") + itoa(count) + _S(")"));

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_FINISH;
  #endif
  }

  TEMPLATE_HEADER
  inline void Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<CONNECT_>)
  {
  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_START("Building of connects : ", 1);
  #endif

    for (unsigned i=0,sz=atomdata_.size(); i<sz; i++)
    {
      unsigned_t &connect_data = atomdata_[i].connect_data;
      connect_data |= 0x00000001; // вставка коннекта с самим собой для действительных атомов
        // полный 0 для коннекта означает dummy атом

    #define INSERT_CONNECTORS(i, matrix) \
      if (matrix.size()) \
      { \
        std::set<int>::iterator it = matrix[i].begin(), ite = matrix[i].end(); \
        for (; it!=ite; ++it) \
        { \
          unsigned n = *it - i; \
          assert(_LT(n, (unsigned)32)); \
          connect_data |= 0x00000001 << n; \
        } \
      }

      INSERT_CONNECTORS(i, t12_matrix_);
      INSERT_CONNECTORS(i, t13_matrix_);
      INSERT_CONNECTORS(i, t14_matrix_);

    #undef INSERT_CONNECTORS
    }

    // Число коннектов должно считаться отдельно, поскольку при вставке
    // коннекта connect_data |= 0x00000001 << n; не проверяется наличие бита
    // в позиции вставки.
    unsigned count = 0;
    for (unsigned i=0,sz=atomdata_.size(); i<sz; i++)
    {
      unsigned_t connect_data = atomdata_[i].connect_data;
      while (connect_data)
      {
        connect_data >>= 1; // сперва удалим коннект с самим собой
        count += connect_data & 0x00000001;
      }
    }

    PRINT_MESSAGE(_S("Have been built -> connects (") + itoa(count) + _S(")"));

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_FINISH;
  #endif
  }

  TEMPLATE_HEADER
  inline void Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<ROTAMER_>)
  {
  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_START("Building of rotamers : ", 1);
  #endif

    unsigned atom_count = atomdata_.size();

    //-----------------------------------------------------------
    //          сделаем массив ротамерных связей
    //-----------------------------------------------------------
    typedef std::pair<unsigned, unsigned>  _Index;

    unsigned sz = bonds_.size();
    std::vector<_Index> bonds__;
    bonds__.reserve(sz);
    for (unsigned i=0; i<sz; i++)
    {
      if (tor12s_.count(i) && !fix12s_.count(i)) continue;
        // исключаем чисто торсионые связи

      const _Bond &bond = bonds_[i];
      unsigned t1 = bond.ndx[0];
      unsigned t2 = bond.ndx[1];
      bonds__.push_back(_Index(t1, t2));
    }

    //-----------------------------------------------------------
    //        сделаем подграфы молекулы (на атомах)
    //-----------------------------------------------------------
    std::vector<int> subgraphs(atom_count);
    unsigned rotamer_count = make_subgraphs(atom_count, &subgraphs[0],
      bonds__.begin(), bonds__.end());

    //-----------------------------------------------------------
    //          сделаем ротамеры молекулы
    //-----------------------------------------------------------
    rotamers_.resize(rotamer_count);
    rotamer_vector_.resize(atom_count); // хранилище номеров атомов для ротамеров
      // 10, 21, 13, 16, 18, ..., 2219 - номера атомов, связанных с ротамерами в сплошном потоке
      // <... rot1 ...> <... rot2 ...> общее количество элементов равно числу атомов

    if (rotamer_count == 1)
    {
      for (unsigned i=0; i<atom_count; i++) rotamer_vector_[i] = i;
      rotamers_[0].set(ATOM, &rotamer_vector_[0]);
    }
    else
    {
      std::vector<unsigned> num_colors(rotamer_count, 0);
      for (unsigned i=0; i<atom_count; i++) num_colors[subgraphs[i]]++;
        // получили число атомов каждого цвета

      std::vector<unsigned> offsets(rotamer_count, 0);
      for (unsigned i=1; i<rotamer_count; i++)
        offsets[i] = offsets[i-1] + num_colors[i-1];

      for (unsigned i=0; i<rotamer_count; i++)
      {
        rotamers_[i].set(ATOM, &rotamer_vector_[offsets[i]]);
        rotamers_[i].count(ATOM) = num_colors[i];
      }

      // запишем атомы в ротамерный массив
      for (unsigned i=0; i<atom_count; i++)
      {
        int color = subgraphs[i]; // цвет атома
        int off = --num_colors[color]; // последний элемент в цвете атома
          // используем число в colors[] как счетчик при обратной записи
          // colors становится недействителен после окончания цикла

        int pos = offsets[color] + off; // позиция ротамера в массиве ротамерных атомов
        rotamer_vector_[pos] = i;
      }

      // сформируем ротамеры
      std::set<unsigned>::iterator sit = tor12s_.begin(), site = tor12s_.end();
      for (; sit != site; ++sit)
      {
        unsigned nb = *sit; // номер парной связи

        if (fix12s_.count(nb)) continue;
          // исключаем невращающиеся торсионые связи

        const _Bond &bond = bonds_[nb];

        unsigned t1 = bond.ndx[0];
        unsigned t2 = bond.ndx[1];
        unsigned c1 = subgraphs[t1];
        unsigned c2 = subgraphs[t2];

        _Rotamer &rotamer1 = rotamers_[c1];
        _Rotamer &rotamer2 = rotamers_[c2];

        unsigned t1__ = 0, t2__ = 0; // вспомогательные атомы, на которые опираются ротамеры
          // при расчете углов между ними, то есть любые из атомов ротамера, кроме стика.

        typedef matrix<ACCUMULATE_>::row_iterator _Iterator;
        _Iterator it = tb_matrix_[t1].begin(), ite = tb_matrix_[t1].end();
        for (; it!=ite; ++it)
        {
          if ((unsigned)*it == nb) continue; // игнорируем связь, совпадающую с ротамером
          const _Bond &bond = bonds_[*it];
          t1__ = bond.ndx[0];
          if (t1__ == t1) t1__ = bond.ndx[1];
          break;
        }

        it = tb_matrix_[t2].begin(), ite = tb_matrix_[t2].end();
        for (; it!=ite; ++it)
        {
          if ((unsigned)*it == nb) continue; // игнорируем связь, совпадающую с ротамером
          const _Bond &bond = bonds_[*it];
          t2__ = bond.ndx[0];
          if (t2__ == t2) t2__ = bond.ndx[1];
          break;
        }

        rotamer1.insert_stick(t2, t1, t1__);
        rotamer2.insert_stick(t1, t2, t2__);
      }
    }

    //----------------------------------------------------------------
    //        сделаем подграфы молекулы (на ротамерах)
    //----------------------------------------------------------------
    bonds__.resize(0); // теперь это связи ротамеров друг с другом

    std::set<unsigned>::iterator sit = tor12s_.begin(), site = tor12s_.end();
    for (; sit != site; ++sit)
    {
      unsigned nb = *sit; // номер парной связи

      if (fix12s_.count(nb)) continue;
        // исключаем невращающиеся торсионые связи

      const _Bond &bond = bonds_[nb];

      unsigned t1 = bond.ndx[0];
      unsigned t2 = bond.ndx[1];
      unsigned r1 = subgraphs[t1];
      unsigned r2 = subgraphs[t2];
      bonds__.push_back(_Index(r1, r2));
    }

    std::vector<int> rotamer_subgraphs(rotamer_count);
    unsigned chain_count = make_subgraphs(rotamer_count, &rotamer_subgraphs[0],
      bonds__.begin(), bonds__.end());
      // число групп связанных ротамеров

    // создадим матрицу, описывающую какие ротамеры попадают в ту или иную цепь
    cr_matrix_.resize(chain_count);
    for (unsigned i=0; i<rotamer_count; i++)
    {
      unsigned color = rotamer_subgraphs[i];
      cr_matrix_[color].insert(i);
    }

    roots_.resize(chain_count, 0);
    std::vector<unsigned> root_amount(chain_count, 0);
    int edge_count = 0; // число ребер во всем (не)связанном графе
    // число ребер графа будет точным при условии, если граф не имеет циклов
    // Циклы могут быть, например, при связывания -S-S- от разных аминокислот

    for (unsigned i=0; i<rotamer_count; i++)
    {
      unsigned color = rotamer_subgraphs[i];
      if (root_amount[color] < rotamers_[i].count(ATOM))
      {
        root_amount[color] = rotamers_[i].count(ATOM);
        roots_[color] = i; // номер ротамера, с которого стартует поиск
      }
      edge_count += rotamers_[i].count(STICK);
    }
    // полное число ребер определяется по формуле  Sum(nsticks) - num_rotamer + nchains
    edge_count -= rotamer_count - chain_count;

    edges_.reserve(edge_count); // резервируем, а не меняем размер,
      // так как точное число ребер неопределено

    //----------------------------------------------------------------
    //       сделаем лес из деревьев по связанным ротамерам
    //----------------------------------------------------------------
    std::vector<bool> used_color(rotamer_count, false); // уже посещенные узлы дерева
    std::stack<unsigned> rotamer_stack; // используем стек как замену рекурсивного прохода

    for (unsigned chain=0; chain<chain_count; chain++)
    {
      rotamer_stack.push(roots_[chain]); // инициализировали стек всеми корнями
      used_color[roots_[chain]] = true;
    }

    while (!rotamer_stack.empty())
    {
      unsigned color1 = rotamer_stack.top();
      rotamer_stack.pop();

      _Rotamer &rotamer1 = rotamers_[color1];
      unsigned count1 = rotamer1.count(STICK);
      for (unsigned stick1=0; stick1<count1; stick1++)
      {
        unsigned n = rotamer1.get(ROTAMER_EXT_ATOM, stick1);
        unsigned color2 = subgraphs[n]; // номер ротамера, к которому принадлежит атом
        if (!used_color[color2]) // избегаем возвратов
        {
          // поиск стика во втором ротамере, соответствующего стику первого ротамера
          _Rotamer &rotamer2 = rotamers_[color2];
          unsigned count2 = rotamer2.count(STICK);
          for (unsigned stick2=0; stick2<count2; stick2++)
          {
            unsigned n = rotamer2.get(ROTAMER_EXT_ATOM, stick2);
            if (subgraphs[n] == (int)color1)
            {
              edges_.push_back(_Edge(color1, stick1, color2, stick2));
              break;
            }
          }

          rotamer_stack.push(color2);
          used_color[color2] = true;
        }
      }
    }
    // Заметим, что массив edges_ согласно вышеприведенному алгоритму является упорядоченным.

  #ifdef FULL_BUILD_DEBUG
    for (unsigned i=0; i<rotamers_.size(); i++) rotamers_[i].print();
  #endif

    PRINT_MESSAGE(_S("Have been built -> rotamers (") + itoa(rotamer_count - roots_.size()) + _S(")"));
      // число связей ротамеров друг с другом (а именно это полагают ротамером) на roots меньше

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_FINISH;
  #endif
  }

  TEMPLATE_HEADER
  INLINE void Archetype_<TEMPLATE_ARG>
  ::build_(_I2T<CHAIN_>)
  {
  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_START("Building of chains : ", 1);
  #endif

    unsigned atom_count = atomdata_.size();

    //-----------------------------------------------------------
    //             сделаем массив связей
    //-----------------------------------------------------------
    typedef std::pair<unsigned, unsigned>  _Index;

    std::vector<_Index> bonds__;
    bonds__.reserve(2 * atom_count);
      // учитываем связи только вперед, что уменьшает запросы на память

    for (unsigned t1=0; t1<atom_count; t1++)
    {
      const _Atomdata &atomdata = atomdata_[t1];
      unsigned nb = atomdata.nbond;
      for (unsigned k=0; k<nb; k++)
      {
        if (atomdata.rnid[k] < 0) continue;
        unsigned t2 = t1 + atomdata.rnid[k];
        bonds__.push_back(_Index(t1, t2));
      }
    }

    //-----------------------------------------------------------
    //        сделаем подграфы молекулы (на атомах)
    //-----------------------------------------------------------
    std::vector<int> subgraphs(atom_count);
    unsigned chain_count = make_subgraphs(atom_count, &subgraphs[0],
      bonds__.begin(), bonds__.end());

    //-----------------------------------------------------------
    //            запишем подграфы молекулы
    //-----------------------------------------------------------
    chains_vector_.resize(atom_count);

    std::vector<unsigned> num_colors(chain_count, 0);
    for (unsigned i=0; i<atom_count; i++) num_colors[subgraphs[i]]++;
      // получили число атомов каждого цвета

    std::vector<unsigned> offsets(chain_count, 0);
    for (unsigned i=1; i<chain_count; i++)
      offsets[i] = offsets[i-1] + num_colors[i-1];

    chains_.resize(chain_count);
    for (unsigned ichain=0; ichain<chain_count; ichain++)
    {
      Chain_ &chain = chains_[ichain];
      chain.set(ATOM, &chains_vector_[offsets[ichain]]);
      chain.count(ATOM) = num_colors[ichain];
    }

    for (unsigned i=0; i<atom_count; i++)
    {
      int color = subgraphs[i]; // цвет атома
      int off = --num_colors[color]; // последний элемент в цвете атома
        // используем число в colors[] как счетчик при обратной записи
        // colors становится недействителен после окончания цикла

      int pos = offsets[color] + off; // позиция ротамера в массиве ротамерных атомов
      chains_vector_[pos] = i;
    }

    PRINT_MESSAGE(_S("Have been built -> chains (") + itoa(chain_count) + _S(")"));

  #ifdef BUILDING_TIME_TESTING
    TIME_TESTING_FINISH;
  #endif
  }

  TEMPLATE_HEADER
  template <typename __Iterator, typename _Iterator>
  inline unsigned Archetype_<TEMPLATE_ARG>
  ::free(_I2T<BOND_>, __Iterator &to, _Iterator it, _Iterator ite) const
  {
    // В реализации используется массив mark для того, чтобы отметить атомы,
    // которые являются свободными. Такая реализация более эффективна, по
    // сравнению с любыми реализация, где используется поиск на соответствие,
    // есть или нет заданный атом в заданной последовательности для данной
    // связи. Заметим, что нельзя использовать внешний массив mark, поскольку
    // это приводит к тому, что функция становится непригодной для параллельного
    // выполнения.

    std::vector<bool> mark(atomdata_.size(), false);
    for (; it!=ite; ++it) mark[*it] = true;

    unsigned count = 0;
    for (unsigned i=0, sz=bonds_.size(); i<sz; i++)
    {
      index_<2> n = bonds_[i].ndx;
      if (mark[n[0]] || mark[n[1]]) { *to++ = i; count++; }
    }

    PRINT_MESSAGE(_S("Have been rebuilt -> free bonds (") + itoa(count) + _S(")"));
    return count;
  }

  TEMPLATE_HEADER
  template <typename __Iterator, typename _Iterator>
  inline unsigned Archetype_<TEMPLATE_ARG>
  ::free(_I2T<ANGLE_>, __Iterator &to, _Iterator it, _Iterator ite) const
  {
    // В реализации используется массив mark для (см. build<FREE_BONDS_>)
    std::vector<bool> mark(atomdata_.size(), false);
    for (; it!=ite; ++it) mark[*it] = true;

    unsigned count = 0;
    for (unsigned i=0, sz=angles_.size(); i<sz; i++)
    {
      index_<3> n = angles_[i].ndx;
      if (mark[n[0]] || mark[n[1]] || mark[n[2]]) { *to++ = i; count++; }
    }

    PRINT_MESSAGE(_S("Have been rebuilt -> free angles (") + itoa(count) + _S(")"));
    return count;
  }

  TEMPLATE_HEADER
  template <typename __Iterator, typename _Iterator>
  inline unsigned Archetype_<TEMPLATE_ARG>
  ::free(_I2T<TORSION_>, __Iterator &to, _Iterator it, _Iterator ite) const
  {
    // В реализации используется массив mark для (см. build<FREE_BONDS_>)
    std::vector<bool> mark(atomdata_.size(), false);
    for (; it!=ite; ++it) mark[*it] = true;

    unsigned count = 0;
    for (unsigned i=0, sz=torsions_.size(); i<sz; i++)
    {
      index_<4> n = torsions_[i].ndx;
      if (mark[n[0]] || mark[n[1]] || mark[n[2]] || mark[n[3]]) { *to++ = i; count++; }
    }

    PRINT_MESSAGE(_S("Have been rebuilt -> free torsions (") + itoa(count) + _S(")"));
    return count;
  }

  TEMPLATE_HEADER
  template <typename __Iterator, typename _Iterator>
  inline unsigned Archetype_<TEMPLATE_ARG>
  ::free(_I2T<PAIR14_>, __Iterator &to, _Iterator it, _Iterator ite) const
  {
    // В реализации используется массив mark для (см. build<FREE_BONDS_>)
    std::vector<bool> mark(atomdata_.size(), false);
    for (; it!=ite; ++it) mark[*it] = true;

    unsigned count = 0;
    for (unsigned i=0, sz=pair14s_.size(); i<sz; i++)
    {
      _Pair14 n = pair14s_[i];
      if (mark[n[0]] || mark[n[1]]) { *to++ = i; count++; }
    }

    PRINT_MESSAGE(_S("Have been rebuilt -> free 1-4 pairs (") + itoa(count) + _S(")"));
    return count;
  }

  TEMPLATE_HEADER
  template <typename _Atom, typename _Iterator>
  inline _E(real_t) Archetype_<TEMPLATE_ARG>
  ::U(_I2T<BOND_>, const _Atom *atoms, _Iterator start, _Iterator end, bool make_print) const
  {
    if (bonds_.size() == 0) return 0.;
    _E(real_t) energy = 0.;

  #ifndef SKIP_BONDS_ENERGY
  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_START("[TIME_TESTING BONDS]", 1)
  #endif

    typedef typename _Bond::param_type   _Param;
    typedef Interaction1_<BOND_>         _Interection;

    _E(real_t) energy__, max_energy = 0.; vector_t R; unsigned count = 0;
    for (; start!=end; ++start, count++)
    {
      const _Bond &bond = bonds_[*start];
      const _Atom &atom0 = atoms[bond.ndx[0]];
      const _Atom &atom1 = atoms[bond.ndx[1]];
      const _Param &param = bond.param;

      R = atom1.X - atom0.X;
      real_t d = (real_t) sqrt(scalar_product(R, R));
      if (d > param.q0 * DEBUG_BOND_FACTOR)
      {
        _S msg = _S("[WARNING] It has been found extra large bond [" ) + make_string(d)
          + _S(" vs ") + make_string(param.q0) + _S("] between atoms: \n")
          + make_string(atom0, atomdata_[bond.ndx[0]]) + _S("\n")
          + make_string(atom1, atomdata_[bond.ndx[1]]) + _S("\n");
        PRINT_MESSAGE(msg);
      }
      if (d < param.q0 / DEBUG_BOND_FACTOR)
      {
        _S msg = _S("[WARNING] It has been found extra small bond [" ) + make_string(d)
          + _S(" vs ") + make_string(param.q0) + _S("] between atoms: \n")
          + make_string(atom0, atomdata_[bond.ndx[0]]) + _S("\n")
          + make_string(atom1, atomdata_[bond.ndx[1]]) + _S("\n");
        PRINT_MESSAGE(msg);
      }
      energy__ = _Interection::U(param.ke, param.q0, R);
      max_energy = fabs(max_energy) < fabs(energy__) ? energy__ : max_energy;
        // замена if на оператор выбора, так как на Cell(?) это быстрее

      energy += energy__;
    }

    if (make_print)
    {
      real_t evarage = (real_t) (count ? energy / count : 0.);
      std::string msg
        = make_string("  U/bond   / (%5d) : %12.5e", count, (float)energy)
        + make_string("   <U> : %12.5e", (float)evarage)
        + make_string("   <maxU> : %12.5e", (float)max_energy);
          // обязательное приведение (real_t)coul_energy требуется, чтобы
          // избежать вывода _E(real_t) типа со спецификатором %12.5le,
          // что неверно в случае real_t == double

      PRINT_MESSAGE(msg);
    }

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_FINISH
  #endif
  #endif // SKIP_BONDS_ENERGY

    return energy;
  }

  TEMPLATE_HEADER
  template <typename _Atom, typename _Iterator>
  inline _E(real_t) Archetype_<TEMPLATE_ARG>
  ::U(_I2T<ANGLE_>, const _Atom *atoms, _Iterator start, _Iterator end, bool make_print) const
  {
    if (angles_.size() == 0) return 0.;
    _E(real_t) energy = 0.;

  #ifndef SKIP_ANGLS_ENERGY
  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_START("[TIME_TESTING ANGLES]", 1)
  #endif

    typedef typename _Angle::param_type  _Param;
    typedef Interaction1_<ANGLE_>        _Interection;

    vector_t A_, B_; _E(real_t) energy__, max_energy = 0.; unsigned count = 0;
    for (; start!=end; ++start, count++)
    {
      const _Angle &angle = angles_[*start];
      const _Atom &atom0 = atoms[angle.ndx[0]];
      const _Atom &atom1 = atoms[angle.ndx[1]];
      const _Atom &atom2 = atoms[angle.ndx[2]];
      const _Param &param = angle.param;

      A_ = atom0.X - atom1.X;
      B_ = atom2.X - atom1.X;
      energy__ = _Interection::U(param.ke, param.q0, A_, B_);

      max_energy = fabs(max_energy) < fabs(energy__) ? energy__ : max_energy;
        // замена if на оператор выбора, так как на Cell(?) это быстрее
      energy += energy__;
    }

    if (make_print)
    {
      real_t evarage = (real_t) (count ? energy / count : 0.);
      std::string msg
        = make_string("  U/angle  / (%5d) : %12.5e", count, (float)energy)
        + make_string("   <U> : %12.5e", (float)evarage)
        + make_string("   <maxU> : %12.5e", (float)max_energy);
      PRINT_MESSAGE(msg);
    }

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_FINISH
  #endif
  #endif

    return energy;
  }

  TEMPLATE_HEADER
  template <typename _Atom, typename _Iterator>
  inline _E(real_t) Archetype_<TEMPLATE_ARG>
  ::U(_I2T<TORSION_>, const _Atom *atoms, _Iterator start, _Iterator end, bool make_print) const
  {
    if (torsions_.size() == 0) return 0.;
    _E(real_t) energy = 0.;

  #ifndef SKIP_TORAS_ENERGY
  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_START("[TIME_TESTING TORSIONS]", 1)
  #endif

    typedef typename _Torsion::param_type  _Param;
    typedef Interaction1_<TORSION_>        _Interection;

    vector_t XA_, XB_, XC_, XD_; _E(real_t) energy__, max_energy = 0.;
    unsigned count = 0;
    for (; start!=end; ++start, count++)
    {
      const _Torsion &torsion = torsions_[*start];
      const _Atom &atom0 = atoms[torsion.ndx[0]];
      const _Atom &atom1 = atoms[torsion.ndx[1]];
      const _Atom &atom2 = atoms[torsion.ndx[2]];
      const _Atom &atom3 = atoms[torsion.ndx[3]];
      const _Param &param = torsion.param;

      energy__ = _Interection::U(param.curf, &param.n[0], &param.v[0], &param.phi[0],
        atom0.X, atom1.X, atom2.X, atom3.X);

      max_energy = fabs(max_energy) < fabs(energy__) ? energy__ : max_energy;
        // замена if на оператор выбора, так как на Cell(?) это быстрее
      energy += energy__;
    }

    if (make_print)
    {
      real_t evarage = (real_t) (count ? energy / count : 0.);
      std::string msg
        = make_string("  U/torsion/ (%5d) : %12.5e", count, (float)energy)
        + make_string("   <U> : %12.5e", (float)evarage)
        + make_string("   <maxU> : %12.5e", (float)max_energy);
      PRINT_MESSAGE(msg);
    }

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_FINISH
  #endif
  #endif

    return energy;
  }

  TEMPLATE_HEADER
  template <typename _Atom, typename _Iterator>
  inline _E(real_t) Archetype_<TEMPLATE_ARG>
  ::U(_I2T<PAIR14_>, const _Atom *atoms, _Iterator start, _Iterator end, bool make_print) const
  {
    if (pair14s_.size() == 0) return 0.;
    _E(real_t) coul_energy = 0., vdw__energy = 0.;

  #if !defined(SKIP_COUL_ENERGY) || !defined(SKIP_VDW_ENERGY)

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_START("[TIME_TESTING PAIRS14]", 1)
  #endif

    real_t coul_scale = forcefield_->get_scale(_I2T<COUL14_>());
    real_t vdw_scale = forcefield_->get_scale(_I2T<VDW14_>());

    vector_t R;
    _E(real_t) max_coul_energy = 0., max_vdw_energy = 0.; unsigned count = 0;
    for (; start!=end; ++start, count++)
    {
      const _Pair14 &pair14 = pair14s_[*start];
      const _Atom &atom0 = atoms[pair14[0]];
      const _Atom &atom1 = atoms[pair14[1]];

      real_t sigma2 = calculate_sigma(atom0.sigma, atom1.sigma);
      real_t eps = atom0.eps * atom1.eps;
      real_t charge = atom0.charge * atom1.charge;
      real_t r2 = distance2(atom0.X, atom1.X);

      if (r2 < sqr(DEBUG_PAIR_FACTOR) * sigma2 || r2 < sqr(MINIMAL_DISTANCE_BETWEEN_ATOMS))
      {
        _S msg = _S("[WARNING] It has been found extra small contact [" )
          + make_string(sqrt(r2)) + _S(" vs ") + make_string(sqrt(sigma2))
          + _S("] between atoms: \n")
          + make_string(atom0, atomdata_[pair14[0]]) + _S("\n")
          + make_string(atom1, atomdata_[pair14[1]]) + _S("\n");
        PRINT_MESSAGE(msg);
      }

      real_t tau2 = r2 / sigma2;
      real_t r_1 = (real_t) (1. / sqrt(r2));

      _E(real_t) coul_energy__ = coul_scale * _Interaction::U(_I2T<COUL_>(), sigma2, charge, r2, tau2, r_1);
      coul_energy += coul_energy__;

      _E(real_t) vdw_energy__  = vdw_scale * _Interaction::U(_I2T<VDW_ >(), sigma2, eps, r2, tau2, r_1);
      vdw__energy += vdw_energy__;

      if (make_print)
      {
        max_coul_energy = fabs(max_coul_energy) < fabs(coul_energy__) ? coul_energy__: max_coul_energy;
        max_vdw_energy = fabs(max_vdw_energy) < fabs(vdw_energy__) ? vdw_energy__: max_vdw_energy;
          // замена if на оператор выбора, так как на Cell(?) это быстрее
      }
    }

    if (make_print)
    {
      #ifndef SKIP_COUL_ENERGY
      {
        real_t coul_evarage = (real_t) (count ? coul_energy / count : 0.);
        std::string msg
          = make_string("  U/coul14 / (%5d) : %12.5e", count, (float)coul_energy)
          + make_string("   <U> : %12.5e", (float)coul_evarage)
          + make_string("   <maxU> : %12.5e", (float)max_coul_energy);
        PRINT_MESSAGE(msg);
      }
      #endif
      #ifndef SKIP_VDW_ENERGY
      {
        real_t vdw_evarage = (real_t) (count ? vdw__energy / count : 0.);
        std::string msg
          = make_string("  U/vdw14  / (%5d) : %12.5e", count, (float)vdw__energy)
          + make_string("   <U> : %12.5e", (float)vdw_evarage)
          + make_string("   <maxU> : %12.5e", (float)max_vdw_energy);
        PRINT_MESSAGE(msg);
      }
      #endif
    }

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_FINISH
  #endif

  #endif // !defined(SKIP_COUL_ENERGY) || !defined(SKIP_VDW_ENERGY)

    return coul_energy + vdw__energy;
  }
/*
  TEMPLATE_HEADER
  template <typename _Atom, typename _Iterator>
  inline _E(real_t) Archetype_<TEMPLATE_ARG>
  ::dU__dQ(_Atom *atoms, _Iterator start, _Iterator end) const
  {
    unsigned int numberi_of_atoms = atomdata_.size();
    for (unsigned i = 0; i < number_of_atoms - 1; i++)
    {
      total_charge += atoms[i].charge;
    }
    atoms[i].charge = -total_
    for (; start != end; start++)
    {
      
    }
  }
*/
  TEMPLATE_HEADER
  template <typename _Atom, typename _Iterator>
  inline _E(real_t) Archetype_<TEMPLATE_ARG>
  ::dU__dX(_I2T<BOND_>, _Atom *atoms, _Iterator start, _Iterator end) const
  {
    if (bonds_.size() == 0) return 0.;
    _E(real_t) energy = 0.;

  #ifndef SKIP_BONDS_ENERGY

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_START("[TIME_TESTING BONDS]", 1)
  #endif

    typedef typename _Bond::param_type  _Param;
    typedef Interaction1_<BOND_>        _Interection;

    vector_t R; _E(real_t) energy__;
    for (; start!=end; ++start)
    {
      const _Bond &bond = bonds_[*start];
      _Atom &atom = atoms[bond.ndx[0]];
      _Atom &atom__ = atoms[bond.ndx[1]];

      R = atom__.X - atom.X;
      energy__ = _Interection::dU__dX(bond.param.ke, bond.param.q0, R);
      energy += energy__;
      atom.F += R;
      atom__.F -= R;
    }

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_FINISH
  #endif
  #endif

    return energy;
  }

  TEMPLATE_HEADER
  template <typename _Atom, typename _Iterator>
  inline _E(real_t) Archetype_<TEMPLATE_ARG>
  ::dU__dX(_I2T<ANGLE_>, _Atom *atoms, _Iterator start, _Iterator end) const
  {
    if (angles_.size() == 0) return 0.;
    _E(real_t) energy = 0.;

  #ifndef SKIP_ANGLS_ENERGY

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_START("[TIME_TESTING ANGLES]", 1)
  #endif

    typedef typename _Angle::param_type  _Param;
    typedef Interaction1_<ANGLE_>        _Interection;

    vector_t A_, B_;
    for (; start!=end; ++start)
    {
      const _Angle &angle = angles_[*start];
      const _Param &param = angle.param;
      _Atom &atom0 = atoms[angle.ndx[0]];
      _Atom &atom1 = atoms[angle.ndx[1]];
      _Atom &atom2 = atoms[angle.ndx[2]];

      A_ = atom0.X - atom1.X;
      B_ = atom2.X - atom1.X;
      energy += _Interection::dU__dX(param.ke, param.q0, A_, B_);
      atom0.F += A_;
      atom1.F -= A_ + B_;
      atom2.F += B_;
    }

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_FINISH
  #endif
  #endif

    return energy;
  }

  TEMPLATE_HEADER
  template <typename _Atom, typename _Iterator>
  inline _E(real_t) Archetype_<TEMPLATE_ARG>
  ::dU__dX(_I2T<TORSION_>, _Atom *atoms, _Iterator start, _Iterator end) const
  {
    if (torsions_.size() == 0) return 0.;
    _E(real_t) energy = 0.;

  #ifndef SKIP_TORAS_ENERGY

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_START("[TIME_TESTING TORSIONS]", 1)
  #endif

    typedef typename _Torsion::param_type  _Param;
    typedef Interaction1_<TORSION_>        _Interection;

    vector_t XA_, XB_, XC_, XD_;
    for (; start!=end; ++start)
    {
      const _Torsion &torsion = torsions_[*start];
      const _Param &param = torsion.param;
      _Atom &atom0 = atoms[torsion.ndx[0]];
      _Atom &atom1 = atoms[torsion.ndx[1]];
      _Atom &atom2 = atoms[torsion.ndx[2]];
      _Atom &atom3 = atoms[torsion.ndx[3]];

      XA_ = atom0.X;
      XB_ = atom1.X;
      XC_ = atom2.X;
      XD_ = atom3.X;
      energy += _Interection::dU__dX(param.curf, &param.n[0], &param.v[0],
        &param.phi[0], XA_, XB_, XC_, XD_);

      atom0.F += XA_;
      atom1.F += XB_;
      atom2.F += XC_;
      atom3.F += XD_;
    }

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_FINISH
  #endif
  #endif

    return energy;
  }

  TEMPLATE_HEADER
  template <typename _Atom, typename _Iterator>
  inline _E(real_t) Archetype_<TEMPLATE_ARG>
  ::dU__dX(_I2T<PAIR14_>, _Atom *atoms, _Iterator start, _Iterator end) const
  {
    if (pair14s_.size() == 0) return 0.;
    _E(real_t) coul_energy = 0., vdw__energy = 0.;

  #if !defined(SKIP_COUL_ENERGY) || !defined(SKIP_VDW_ENERGY)

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_START("[TIME_TESTING PAIRS14]", 1)
  #endif

    real_t coul_scale = forcefield_->get_scale(_I2T<COUL14_>());
    real_t vdw_scale = forcefield_->get_scale(_I2T<VDW14_>());

    _E(real_t) coul_energy__, vdw__energy__;
    vector_t R;

    for (; start!=end; ++start)
    {
      const _Pair14 &pair14 = pair14s_[*start];
      _Atom &atom0 = atoms[pair14[0]];
      _Atom &atom1 = atoms[pair14[1]];

      real_t sigma2 = calculate_sigma(atom0.sigma, atom1.sigma);
      real_t eps = atom0.eps * atom1.eps;
      real_t charge = atom0.charge * atom1.charge;

      R = atom1.X - atom0.X;
      real_t r2 = scalar_product(R, R);

      real_t tau2 = r2 / sigma2;
      real_t r_1 = (real_t) (1. / sqrt(r2));

      real_t coef_du__dq = (real_t) 0.;

      real_t coul_du__dq = (real_t) 0.;
      coul_energy__ = _Interaction::dU__dX(_I2T<COUL_>(), &coul_du__dq, sigma2, charge, r2, tau2, r_1);
      coul_energy += coul_scale * coul_energy__;
      coef_du__dq += coul_scale * coul_du__dq;

      real_t vdw__du__dq = (real_t) 0.;
      vdw__energy__ = _Interaction::dU__dX(_I2T<VDW_ >(), &vdw__du__dq, sigma2, eps, r2, tau2, r_1);
      vdw__energy += vdw_scale * vdw__energy__;
      coef_du__dq += vdw_scale * vdw__du__dq;

      R *= coef_du__dq * r_1;
      atom0.F += R;
      atom1.F -= R;
    }

  #ifdef ENERGY_TIME_TESTING
    TIME_TESTING_FINISH
  #endif

  #endif // !defined(SKIP_COUL_ENERGY) || !defined(SKIP_VDW_ENERGY)

    return coul_energy + vdw__energy;
  }

  #undef TEMPLATE_HEADER
  #undef TEMPLATE_ARG
}
#endif
