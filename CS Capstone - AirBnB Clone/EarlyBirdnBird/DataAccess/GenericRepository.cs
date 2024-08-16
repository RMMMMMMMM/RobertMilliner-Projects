using Infrastructure.Interfaces;
using Infrastructure.Models;
using Microsoft.EntityFrameworkCore;
using System.Linq.Expressions;



namespace DataAccess
{
    //reminder inheritance requires all methods to be included
    public class GenericRepository<T> : IGenericRepository<T> where T : class
    {
        //injects repository we are using communicaiton for DB context which communicates with physical DB
        //using read only appDbContest we cann access to only querying cababilities of DBContext. UnitofWork actuall actually 
        //(Commits to the PHYSICAL Tables (not internal object). 

        private readonly ApplicationDbContext _dbContext;


        public GenericRepository(ApplicationDbContext dbContext)
        {
            _dbContext = dbContext;
        }


        public void Add(T entrity)
        {
            //using framwork to do the add method. the framework converts to the applicable db language through DB context
            _dbContext.Set<T>().Add(entrity);
            
        }

        public void Delete(T entrity)
        {
            _dbContext.Set<T>().Remove(entrity);

        }

        public void Delete(IEnumerable<T> entrity)
        {   
            //use remove range for the Ienumerable list
            _dbContext.Set <T>().RemoveRange(entrity);

        }

        //virtual keyword is used to modify the method property indxer or
        //allows for it to be overridden in a derived class
        public virtual T Get(Expression<Func<T, bool>> predicate, bool trackChanges = false, string? includes = null)
        {
            //determine if we are joing other objects
            if (includes == null)
            {
                if (!trackChanges) //if fallse
                {
                    return _dbContext.Set<T>().Where(predicate).AsNoTracking().FirstOrDefault();
                }
                else //we tacking changes entity freamwork does by default
                {
                    return _dbContext.Set<T>().Where(predicate).FirstOrDefault();
                }
            }
            else // have includes
            {
                //includes =comma seperated objs without spaces

                //cast data set a s"Iqueryable" so that we can iterate line by line
                IQueryable<T> queryable = _dbContext.Set<T>();

                //spits includes to new char array when "," is seen
                foreach (var includeProperty in includes.Split(new char[] {','}, StringSplitOptions.RemoveEmptyEntries))
                {
                    queryable = queryable.Include(includeProperty); //thiis joins tables based on primary key

                }

                if (!trackChanges) //if fallse
                {
                    return queryable.Where(predicate).AsNoTracking().FirstOrDefault();
                }
                else //we tacking changes entity freamwork does by default
                {
                    return queryable.Where(predicate).FirstOrDefault();
                }
                
            }
        }

        public IEnumerable<T> GetAll(Expression<Func<T, bool>>? predicate = null, Expression<Func<T, int>>? orderBy = null, string? includes = null)
        {
            //cast data set a s"Iqueryable" so that we can iterate line by line
            IQueryable<T> queryable = _dbContext.Set<T>();


            //determine if we are joing other objects
            if (predicate != null && includes == null)
            {
                return _dbContext.Set<T>().Where(predicate).AsEnumerable();
            }
            else if (includes != null)// have includes
            {
                //spits includes to new char array when "," is seen
                foreach (var includeProperty in includes.Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries))
                {
                    queryable = queryable.Include(includeProperty); //thiis joins tables based on primary key
                }
            }

            if(predicate == null)
            {
                if (orderBy == null)
                {
                    return queryable.AsEnumerable();
                }
                else
                {
                    return queryable.OrderBy(orderBy).ToList();
                }
            }
            else
            {
                if(orderBy == null)
                {
                    return queryable.Where(predicate).ToList(); 
                }
                else
                {
                    return queryable.Where(predicate).OrderBy(orderBy).ToList();
                }
            }

        }

        public async Task<IEnumerable<T>> GetAllAsync(Expression<Func<T, bool>>? predicate = null, Expression<Func<T, int>>? orderBy = null, string? includes = null)
        {
            //cast data set a s"Iqueryable" so that we can iterate line by line
            IQueryable<T> queryable = _dbContext.Set<T>();

            //determine if we are joing other objects
            if (predicate != null && includes == null)
            {
                return _dbContext.Set<T>().Where(predicate).AsEnumerable();
            }
            else if (includes != null)// have includes
            {
                //spits includes to new char array when "," is seen
                foreach (var includeProperty in includes.Split(new char[] { ',' }, StringSplitOptions.RemoveEmptyEntries))
                {
                    queryable = queryable.Include(includeProperty); //thiis joins tables based on primary key
                }
            }

            if (predicate == null)
            {
                if (orderBy == null)
                {
                    return queryable.AsEnumerable();
                }
                else
                {
                    return await queryable.OrderBy(orderBy).ToListAsync();
                }
            }
            else
            {
                if (orderBy == null)
                {
                    return await queryable.Where(predicate).ToListAsync();
                }
                else
                {
                    return await queryable.Where(predicate).OrderBy(orderBy).ToListAsync();
                }
            }
        }

        //add await any time we return async
        public virtual async Task<T> GetAsync(Expression<Func<T, bool>> predicate, bool trackChanges = false, string? includes = null)
        {
            //determine if we are joing other objects
            if (includes == null)
            {
                if (!trackChanges) //if fallse
                {
                    return await _dbContext.Set<T>().Where(predicate).AsNoTracking().FirstOrDefaultAsync();
                }
                else //we tacking changes entity freamwork does by default
                {
                    return await _dbContext.Set<T>().Where(predicate).FirstOrDefaultAsync();
                }
            }
            else // have includes
            {
                //includes =comma seperated objs without spaces

                //cast data set a s"Iqueryable" so that we can iterate line by line
                IQueryable<T> queryable = _dbContext.Set<T>();

                //spits includes to new char array when "," is seen
                foreach (var includeProperty in includes.Split(new char[] {','}, StringSplitOptions.RemoveEmptyEntries))
                {
                    queryable = queryable.Include(includeProperty); //thiis joins tables based on primary key

                }

                if (!trackChanges) //if fallse
                {
                    return queryable.Where(predicate).AsNoTracking().FirstOrDefault();
                }
                else //we tacking changes entity freamwork does by default
                {
                    return queryable.Where(predicate).FirstOrDefault();
                }

            }
        }


        //virtual keyword is used to modify the method property indxer or
        //allows for it to be overridden in a derived class
        public virtual T GetById(int? id)
        {
            return _dbContext.Set<T>().Find(id);
        }

        public void Update(T entrity)
        {
            //tracks changes, flagges record as modified to the system
            _dbContext.Entry(entrity).State = EntityState.Modified;

        }



    }
}
